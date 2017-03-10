Snap is a cross-distribution packaging format for Linux. See
https://snapcraft.io for more information.

How to build
------------
Just type
```
snapcraft
```
from the *parent* directory: that is, not from the directory containing this
README file, but from its parent directory. This snap has been successfully
built in Ubuntu 16.04; more recent Ubuntu versions should work as well. Note
that, regardless of the environment used to build it, the generated snap
package will work on any distribution where snap is installed.

How to test the snap
--------------------
The generated snap can be installed by typing
```
snap install --dangerous ./openmvg*.snap
```
The `--dangerous` flag is needed because the snap has not been verified by the
store.

How to upload the snap to the store
-----------------------------------
A thorough guide can be found at the [snapcraft.io
site](https://snapcraft.io/docs/build-snaps/publish).

