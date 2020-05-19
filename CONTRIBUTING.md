# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, slack, or any other method with the owners of this repository before making a change. 

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Main Process

- Follow [GitHub flow](https://guides.github.com/introduction/flow/)
- [Commit with emoji](https://gist.github.com/parmentf/035de27d6ed1dce0b36a)
- Follow repo's branch naming convention:
  - [Feature Branches](https://gist.github.com/digitaljhelms/4287848#feature-branches):
    - `features/feature-name`
    - `features/feature-area/feature-name`
  - [Bug branches](https://gist.github.com/digitaljhelms/4287848#bug-branches): `bugfix/description`
  - [Hotfix branches](https://gist.github.com/digitaljhelms/4287848#hotfix-branches): `hotfix/description`
- Use repo's [Pull Request Template](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/master/.github/PULL_REQUEST_TEMPLATE.md)
- Rebase merge is required avoid commit too many updates histories
- Require two reviewers' approval for a merge
- Requester should merge their own Pull Request once it gets sign-off unless they don't have permission to do so


## Code Review Rule
1. Each PR review must have at least 2 same-team developer approvals. Manager approval does not count.
2. Each PR must have a good description. From reading the description, the reviewer should be able to understand what the code is meant to do. This has to be true even if there is a ticket or a requirements page.
3. PR must have sufficient unit test and integration test coverage.
4. If the PR is a bug fix, it must contain a test such that, should the bug fix be reverted, this test would fail.

## Code Style

### Python
[`Pep-8`](https://www.python.org/dev/peps/pep-0008/#tabs-or-spaces) should be
followed for code formatting at a minimum.
The [`black`](https://github.com/psf/black) formatter should be used for any
new repositories and its use is suggested for formatting any new changes or
additions to repositories that do not yet have project-wide conformance with
black.
Most repositories will and should check for formatting using
[`pycodestyle`](https://pypi.org/project/pycodestyle/) as part of the CI checks.
```
black {source_file_or_directory}
```


### R
[Google R style guide](https://google.github.io/styleguide/Rguide.html) should be
followed for code formatting at a minimum. It is not important to squabble about what is a good style. It is absolutely essential to have a style, and then to be 100% consistent with it.
Code should be checked by [`lintr`](https://www.rdocumentation.org/packages/lintr/versions/1.0.3) before committing.
```R
lintr::lint("/path/to/Rscript.R")
```

### CWL
[D3B CWL style guide](https://github.com/d3b-center/handbook/blob/DEV-WIP/docs/development/cwl_style_guide.md) should be followed for any CWL projects.


### Bash/Shell
We recommend staying away from shell scripts as much as possible. A language like R or Python is almost always a better choice. Use shell scripts only if there’s a strong restriction on project’s dependencies size or any other requirements that are more important in a particular case.

[`ShellCheck`](https://www.shellcheck.net/) should be used to lint our shell scripts.
```shell
shellcheck {source_file}
```

[Google Shell Style Guide](https://google.github.io/styleguide/shell.xml) should be followed, and [`shfmt`](https://formulae.brew.sh/formula/shfmt) is recommended to maintain consistent formatting. 
```shell
shfmt -i 2 -ci -w scripts/**/*.sh

-i uint   indent: 0 for tabs (default), >0 for number of spaces
-ci       switch cases will be indented
-w        write result to file instead of stdout
```


## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at [INSERT EMAIL ADDRESS]. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
