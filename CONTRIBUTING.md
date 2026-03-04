# How to Contribute

`ClaMS: Clustering Comparison` (`ClaMS-CC` throughout) is an open source project 
that is a stand-alone part of the [ClaMS](https://github.com/llnl/clams) project.
Our team welcomes contributions from collaborators in the form of raising issues
as well as code contributions including hotfixes, code improvements, and new
features.

`ClaMS-CC` is distributed under the terms of the
[BSD-Commercial license](https://github.com/LLNL/clams-cc/blob/main/LICENSE).
All new contributions must be made under this license.

If you identify a problem such as a bug or awkward or confusing code, or require
a new feature, please feel free to start a thread on our
[issue tracker](https://github.com/LLNL/clams-cc/issues).
Please first review the existing issues prior to avoid duplicate issues.

If you plan on contributing to `ClaMS-CC`, please review the
[issue tracker](https://github.com/LLNL/clams-cc/issues) to check for
threads related to your desired contribution.
We recommend creating an issue prior to issuing a pull request if you are
planning significant code changes or have questions.

# Contribution Workflow

These guidelines assume that the reader is familiar with the basics of
collaborative development using git and GitHub.
This section will walk through our preferred pull request workflow for
contributing code to `ClaMS-CC`.
While releases are tagged `v0.XX` (e.g., `v0.3`) on the main branch, with larger
numbers indicating more recent releases, the current main development branch
will always be titled `v0.XX-dev` (.e.g., `v0.4-dev`).
While we will periodically delete old dev branches, if there are multiple
`v0.XX-dev` branches then the larger numbered branch will be the current
development branch.
The tl;dr guidance is:
- Fork the [LLNL ClaMS-CC repository](https://github.com/LLNL/clams-cc)
- Create a descriptively named branch
(`feature/myfeature`, `iss/##`, `hotfix/bugname`, etc) in your fork off of
the most current `v0.XX-dev` branch
- Commit code, following our [guidelines](#formatting-guidelines)
- Create a [pull request](https://github.com/LLNL/clams-cc/compare) from
your branch targeting the most current LLNL `v0.XX-dev` branch

## Forking ClaMS-CC

If you are not a `ClaMS` developer at LLNL, you will not have permissions to push
new branches to the repository.
Even `ClaMS` developers at LLNL will want to use forks for most contributions.
This will create a clean copy of the repository that you own, and will allow for
exploration and experimentation without muddying the history of the central
repository.

If you intend to maintain a persistent fork of `ClaMS-CC`, it is a best practice to
set the LLNL repository as the `upstream` remote in your fork.
```
$ git clone git@github.com:your_name/clams-cc.git
$ cd clams-cc
$ git remote add upstream git@github.com:LLNL/clams-cc.git
```
This will allow you to incorporate changes to the `main` and `v0.XX-dev`
branches as they evolve.
For example, to your fork's `v0.XX-dev` branch perform the following commands:
```
$ git fetch upstream
$ git checkout v0.XX-dev
$ git pull upstream v0.XX-dev
$ git push origin v0.XX-dev
```
It is important to keep your develop branch up-to-date to reduce merge conflicts
resulting from future PRs.

## Contribution Types

Most contributions will fit into one of the follow categories, which by
convention should be committed to branches with descriptive names.
Here are some examples:
- A new feature (`feature/<feature-name>`)
- A bug or hotfix (`hotfix/<bug-name>` or `hotfix/<issue-number>`)
- A response to a [tracked issue](https://github.com/LLNL/clams-cc/issues)
(`iss/<issue-number>`)
- A work in progress, not to be merged for some time (`wip/<change-name>`)

### Developing a new feature

New features should be based on the `v0.XX-dev` branch:
```
$ git checkout v0.XX-dev
$ git pull upstream v0.XX-dev
```
You can then create new local and remote branches on which to develop your
feature.
```
$ git checkout -b feature/<feature-name>
$ git push --set-upstream origin feature/<feature-name>
```
Commit code changes to this branch, and add any tests to `${CLAMS_CC_ROOT}/tests`
that validate the correctness of your code, modifying existing tests if need be.
Be sure to add `add_clams_cc_mpi_test(<test-name>)` or `add_clams_cc_test(<test-name>)` 
(if MPI is not needed) to `${CLAMS_CC_ROOT}/tests/CMakeLists.txt` as appropriate.
Be sure that you can build and that `make test` passes, and that your test
runs successfully.
Once testing is implemented, branches with PRs targeting `v0.XX-dev` branches
will trigger CI jobs upon pushes.

Make sure that you follow our [formatting guidelines](#formatting-guidelines)
for any changes to the source code or build system.
If you create new methods or classes, please add Doxygen documentation
(guidelines forthcoming).

Once your feature is complete and your tests are passing, ensure that your
remote fork is up-to-date and
[create a PR](https://github.com/LLNL/clams-cc/compare).

### Developing a hotfix

Firstly, please check to ensure that the bug you have found has not already been
fixed in `v0.XX-dev`.
If it has, we suggest that you temporarily swap to the `v0.XX-dev` branch until
the next release.

If you have identified an unsolved bug, you can document the problem and create
an [issue](https://github.com/LLNL/clams-cc/issues).
If you would like to solve the bug yourself, follow a similar protocol to
feature development.
First, ensure that your fork's `v0.XX-dev` branch is up-to-date.
```
$ git checkout v0.XX-dev
$ git pull upstream v0.XX-dev
```
You can then create new local and remote branches on which to write your bug
fix.
```
$ git checkout -b hotfix/<bug-name>
$ git push --set-upstream origin hotfix/<bug-name>
```

Firstly, create a test added to `${CLAMS_CC_ROOT}/tests` that reproduces the bug or
modify an existing test to catch the bug if that is more appropriate.
Be sure to add `add_clams_cc_mpi_test(<test-name>)` or `add_clams_cc_test(<test-name>)`
to `${CLAMS_CC_ROOT}/tests/CMakeLists.txt` if you create a new test case.
Then, modify the code to fix the bug and ensure that your new or modified test
case(s) pass via `make test`.

Please update function and class documentation to reflect any changes as
appropriate, and follow our [formatting guidelines](#formatting-guidelines) with
any new code.

Once you are satisfied that the bug is fixed, ensure that your remote fork is
up-to-date and [create a PR](https://github.com/LLNL/clams-cc/compare).

# Tests

While this is not currently implemented, `ClaMS-CC` will use GitHub actions for
continuous integration tests.
Our tests will run automatically against every new commit and pull request on
PRs targeting a dev branch, and pull requests must pass all tests prior to being 
considered for merging into the main project. If you are developing a new feature 
or fixing a bug, please add a test or modify existing tests that will ensure the 
correctness of the new code.

`ClaMS-CC`'s tests are contained in the `test` directory, and new tests must be added
manually to `test/CMakeLists.txt` by adding `add_clams_cc_mpi_test(<test-name>)` or
`add_clams_cc_test(<test-name>)` in order to ensure that they are automatically
built and checked via `make test`.

# Formatting Guidelines

## Code Style

`ClaMS-CC` uses
[clang-format](https://www.kernel.org/doc/html/v4.17/process/clang-format.html)
to guarantee a consistent format for C++ code.
Our style settings are located in
[.clang-format](https://github.com/LLNL/clams-cc/blob/master/.clang-format)
in the project root.
`clang-format` is easy to use, and can be easily instrumented to auto-format
code using most modern editors:
- [vscode](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format)
- [sublime](https://packagecontrol.io/packages/Clang%20Format)
- [vi](https://github.com/rhysd/vim-clang-format)
- [emacs](https://github.com/sonatard/clang-format)

`ClaMS-CC` also uses [cmake-format](https://github.com/cheshirekow/cmake_format) to
guarantee a consistent format for cmake build code.
Our style settings are located in .cmake-format.py in the project root.