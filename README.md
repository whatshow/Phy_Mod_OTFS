# Mod_OTFS
OTFS modulation

* **In another local repositiory, add this module**
```sh
git submodule add git@github.com:whatshow/Phy_Mod_OTFS.git Modules/Mod_OTFS
```
Now, you can see a folder `Modules` with `Mod_OTFS` inside

* **Import this model**
    * Matlab
        ```
        addpath("Modules/OTFS");
        ```
    * Python
        ```python
        if '.' not in __name__ :
            from Modules.Mod_OTFS.OTFS import OTFS
        else:
            from .Modules.Mod_OTFS.OTFS import OTFS
        ```

## How to use
