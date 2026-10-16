# IAM for the image workflows

`build-image.yml` and `promote-image.yml` publish to the `swx-trec` public
registry (`public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa`), which belongs to the
dev account (174369828756). Each assumes its own role there through GitHub's
OIDC provider, so this repository holds no AWS keys. The EC2 instances pull
the image without credentials, so nothing is needed in any other account.

These documents are the roles, kept here so they can be rebuilt and reviewed.

| role | assumed by | can |
|---|---|---|
| `GitHubH3liovizImagePublisher` | `build-image.yml` | log in to ECR Public, push images to `pvw-h3lioviz-osmesa`, read their digests |
| `GitHubH3liovizImagePromoter` | `promote-image.yml` | read digests in `pvw-h3lioviz-osmesa` and point tags at images already there |

Neither role can delete an image. Both can call `PutImage`, which can point
any tag at another image; the workflows only ever write `sha-*` and `dev`
(build) or `prod` and `prod-until-*` (promote), never `v1`, `noaa` or
`noaatest`. That limit is in the workflow files, not in IAM. An image whose
last tag moves away is still in the repository and can be pulled by digest.

## Creating them

Both roles were created in the dev account with these commands on 2026-10-16.
To rebuild them, run from the repository root, with credentials that can write
IAM in the dev account:

```bash
aws iam create-role \
  --role-name GitHubH3liovizImagePublisher \
  --description "Publishes the h3lioviz server image from GitHub Actions" \
  --assume-role-policy-document file://.github/iam/github-h3lioviz-image-publisher.trust.json
aws iam put-role-policy \
  --role-name GitHubH3liovizImagePublisher \
  --policy-name PublishH3liovizServerImage \
  --policy-document file://.github/iam/github-h3lioviz-image-publisher.policy.json

aws iam create-role \
  --role-name GitHubH3liovizImagePromoter \
  --description "Moves the h3lioviz server image prod tag from GitHub Actions" \
  --assume-role-policy-document file://.github/iam/github-h3lioviz-image-promoter.trust.json
aws iam put-role-policy \
  --role-name GitHubH3liovizImagePromoter \
  --policy-name PromoteH3liovizServerImage \
  --policy-document file://.github/iam/github-h3lioviz-image-promoter.policy.json
```

The OIDC provider already exists in the dev account
(`arn:aws:iam::174369828756:oidc-provider/token.actions.githubusercontent.com`,
also used by liveswx_dst and swxtrec-cdk).

## Who can assume them

Only the named workflow file, running from `refs/heads/main` of
SWxTREC/h3lioviz-server. Every claim is matched with `StringEquals`:

- `sub` for the publisher is `repo:SWxTREC/h3lioviz-server:ref:refs/heads/main`.
  Pull request runs and other branches get a different `sub` and are refused,
  which is why the build workflow publishes only on main.
- `sub` for the promoter is `repo:SWxTREC/h3lioviz-server:environment:prod`,
  because a job that names an environment gets that form ("The subject claim
  includes the branch name of the workflow, but only if the job doesn't
  reference an environment", GitHub's OIDC reference). The ref is pinned by
  `job_workflow_ref` instead, and by the environment's branch rule below.
- `job_workflow_ref` pins the workflow file and the branch it runs from.
  Without it any workflow on main that asks for an OIDC token could assume
  the role, since `sub` is the same for every workflow on one ref. (GitHub
  documents this claim for reusable workflows; for a top-level workflow it
  holds the workflow's own path, which liveswx_dst's roles rely on in
  practice.)
- `repository_id` (368988856) and `repository_owner_id` (60273811, the SWxTREC
  organization) pin the repository by GitHub's immutable ids. The names in
  `sub` and `job_workflow_ref` could be taken over if the repository or
  organization were renamed or deleted; AWS recommends the ids for that reason
  ("Use immutable identifiers, not names", IAM condition keys reference).

Renaming or moving either workflow file breaks its authentication until the
trust document here is updated and applied.

## The prod environment (approval before each promotion)

`promote-image.yml` runs in the `prod` environment, so every promotion waits
for a reviewer. It was created on 2026-10-16 with reviewers bmcclellan-cu and
greglucas (either one can approve, including a promotion they started
themselves), administrators unable to bypass it, and `main` as the only
allowed branch. To rebuild it, a repository admin uses Settings, Environments:

1. New environment, named `prod`.
2. Required reviewers: add the people or teams who may approve a promotion
   (up to 6; one approval is enough). Optionally tick "Prevent self-review",
   so whoever starts a promotion cannot approve it.
3. Turn off the option that lets administrators bypass the protection rules.
   It is on by default ("By default, administrators can bypass the protection
   rules and force deployments to specific environments", GitHub's
   environments reference); the API reports it as `can_admins_bypass`.
4. Deployment branches and tags: "Selected branches and tags", with one rule
   for the branch `main`.

The workflow's first step reads the environment back and stops, before any
AWS credentials, unless it has required reviewers, does not let
administrators bypass them, and allows the branch `main` only. That matters
because GitHub creates a missing environment, with no protection, as soon as
a workflow names it.

### The gate is only as strong as `main`

Both roles trust whatever the workflow files on `main` say. `main` requires
a pull request with one approving review, and new commits dismiss earlier
approvals, so a changed workflow needs a second person before it runs with
either role. Administrators are exempt (the protection does not apply to
admins): an administrator can push a changed workflow straight to `main`, and
an edited `build-image.yml` then runs with the publisher role, whose
`PutImage` can point any tag (`prod` included) at any image, with no
environment and no reviewer. Applying the protection to administrators too
would close that.

## Changing them

Edit the JSON here, then apply it, for example:

```bash
aws iam update-assume-role-policy --role-name GitHubH3liovizImagePublisher \
  --policy-document file://.github/iam/github-h3lioviz-image-publisher.trust.json
aws iam put-role-policy --role-name GitHubH3liovizImagePublisher \
  --policy-name PublishH3liovizServerImage \
  --policy-document file://.github/iam/github-h3lioviz-image-publisher.policy.json
```

## Old images

Older images are kept for rollback and for checking what ran when. Each build
gets a tag no other build reuses (`sha-<12 hex of the commit>`, or
`sha-<12 hex>-run<run id>.<attempt>` when the same commit is built again),
and each promotion that moves `prod` gives the image it named before a
`prod-until-<UTC time>` tag. Nothing deletes images automatically.

Very old images can be removed by hand when they are no longer useful, by
someone with delete rights in the dev account. Never remove an image that a
`dev`, `prod`, `v1`, `noaa` or `noaatest` tag names, or the newest
`prod-until-` images, which are the rollback targets.
