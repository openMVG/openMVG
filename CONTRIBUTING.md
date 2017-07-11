# Contributing to OpenMVG

First off all, thanks for taking the time to contribute!

**Contributions** (bug reports or bug fixes, features, etc.) and **feedback**  are welcome and should be submitted thanks to issues form and/or pull requests on GitHub.



Here the list of things to check for contribution:

1. Start your dev from **develop** branch:

  -  Master is storing the current release.
  - Develop store the WIP for the next release candidate. Please submit any new features to the develop branch.

2. Consider adding unit test:
 
  - Add unit tests for all newly added code.
  - It helps to check accuracy of your code implementation and show how to use the API. It helps also any future maintenance.

2. Consider adding documentation:
 
  - Please use doxygen documentation.

     - You can also consider extending the exiting Sphinx documentation (see ./docs/sphinx)


4. Push your change and make a PR towards develop:

  - please check the status of the continuous integration tools:
      -  Compilation (travis, appveyor),
      -  Code quality and style (codacity)).
      
If you want to learn more about OpenSource we refer you to this:
https://opensource.guide/
