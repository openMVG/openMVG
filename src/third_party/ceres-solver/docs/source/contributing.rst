.. _chapter-contributing:

============
Contributing
============

We welcome contributions to Ceres, whether they are new features, bug
fixes or tests. The Ceres `mailing
<http://groups.google.com/group/ceres-solver>`_ list is the best place
for all development related discussions. Please consider joining
it. If you have ideas on how you would like to contribute to Ceres, it
is a good idea to let us know on the mailing list before you start
development. We may have suggestions that will save effort when trying
to merge your work into the main branch. If you are looking for ideas,
please let us know about your interest and skills and we will be happy
to make a suggestion or three.

We follow Google's `C++ Style Guide
<http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml>`_ and
use `git <http://git-scm.com/>`_ for version control. We use the
`Gerrit <https://ceres-solver-review.googlesource.com/>`_ to collaborate and
review changes to Ceres. Gerrit enables pre-commit reviews so that
Ceres can maintain a linear history with clean, reviewed commits, and
no merges.

We now describe how to set up your development environment and submit
a change list for review via Gerrit.

Setting up your Environment
===========================

1. Download and configure ``git``.

   * Mac ``brew install git``.
   * Linux ``sudo apt-get install git``.
   * Windows. Download `msysgit
     <https://code.google.com/p/msysgit/>`_, which includes a minimal
     `Cygwin <http://www.cygwin.com/>`_ install.

2. Sign up for `Gerrit
   <https://ceres-solver-review.googlesource.com/>`_. You will also
   need to sign the Contributor License Agreement (CLA) with Google,
   which gives Google a royalty-free unlimited license to use your
   contributions. You retain copyright.

3. Clone the Ceres Solver ``git`` repository from Gerrit.

   .. code-block:: bash

      git clone https://ceres-solver.googlesource.com/ceres-solver


4. Build Ceres, following the instructions in
   :ref:`chapter-installation`.

   On Mac and Linux, the ``CMake`` build will download and enable
   the Gerrit pre-commit hook automatically. This pre-submit hook
   creates `Change-Id: ...` lines in your commits.

   If this does not work OR you are on Windows, execute the
   following in the root directory of the local ``git`` repository:

   .. code-block:: bash

      curl -o .git/hooks/commit-msg https://ceres-solver-review.googlesource.com/tools/hooks/commit-msg
      chmod +x .git/hooks/commit-msg

5. Configure your Gerrit password with a ``.gitcookies`` which allows pushing
   to Gerrit without having to enter a very long random password every time:

   * Sign into `http://ceres-solver-review.googlesource.com
     <http://ceres-solver-review.googlesource.com>`_.

   * Click ``Settings -> HTTP Password -> Obtain Password``.

   * (maybe) Select an account for multi-login. This should be the
     same as your Gerrit login.

   * Click ``Allow access`` when the page requests access to your
     ``git`` repositories.

   * Follow the instructions from Gerrit to create a ``.gitcookies`` file on
     your system, either in ``$HOME/.gitcookies`` (Mac and Linux) or
     ``%USERPROFILE%\.gitcookies`` (Windows). Note that for Windows, please get
     a recent `Git for Windows <https://git-scm.com/download/win>`_ install to
     enable automatic lookup in the ``%USERPROFILE%\.gitcookies``.

Submitting a change
===================

1. Make your changes against master or whatever branch you
   like. Commit your changes as one patch. When you commit, the Gerrit
   hook will add a `Change-Id:` line as the last line of the commit.

   Make sure that your commit message is formatted in the `50/72 style
   <http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html>`_.

2. Push your changes to the Ceres Gerrit instance:

   .. code-block:: bash

      git push origin HEAD:refs/for/master

   When the push succeeds, the console will display a URL showing the
   address of the review. Go to the URL and add at least one of the
   maintainers (Sameer Agarwal, Keir Mierle, Alex Stewart or William
   Rucklidge) as reviewers.

3. Wait for a review.

4. Once review comments come in, address them. Please reply to each
   comment in Gerrit, which makes the re-review process easier. After
   modifying the code in your ``git`` instance, *don't make a new
   commit*. Instead, update the last commit using a command like the
   following:

   .. code-block:: bash

      git commit --amend -a

   This will update the last commit, so that it has both the original
   patch and your updates as a single commit. You will have a chance
   to edit the commit message as well. Push the new commit to Gerrit
   as before.

   Gerrit will use the ``Change-Id:`` to match the previous commit
   with the new one. The review interface retains your original patch,
   but also shows the new patch.

   Publish your responses to the comments, and wait for a new round
   of reviews.
