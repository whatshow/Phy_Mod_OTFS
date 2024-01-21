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
All OTFS codes are uniform in matlab and python as a class of `OTFS`. This class is the whole process of OTFS. This section will illustrate the methods of this class following the process from Tx to Rx.
* OTFS<br>
    `@nSubcarNum`: the subcarrier number<br>
    `@nTimeslotNum`: the timeslot number<br>
    `@batch_size`**(optional)** : the batch size **(only used in python)**.<br>
    ```matlab
    % matlab
    nSubcarNum = 16;
    nTimeslotNum = 7;
    otfs = OTFS(nSubcarNum, nTimeslotNum);
    ```
    ```python
    # python
    nSubcarNum = 16;
    nTimeslotNum = 7;
    otfs = OTFS(nSubcarNum, nTimeslotNum);
    otfs = OTFS(nSubcarNum, nTimeslotNum, batch_size=4); # using batch
    ```
* modulate<br>
    `@symbols`: a vector of [(batch_size), nSubcarNum*nTimeslotNum], or a matrix in delay Doppler domain [(batch_size), Doppler, delay] or [(batch_size), nTimeslotNum ,nSubcarNum]<br>
    ```matlab
    % matlab
    otfs.modulate(x_origin);
    ```
    ```python
    # python
    otfs.modulate(x_origin);
    ```
* setChannel<br>
    This methods has 2 kinds of inputs.
    * generate random channel<br>
        If we use `batch`, the channel will be different for each sample in the batch.
        `@p`: the path number (scalar)<br>
        `@lmax`: the maximal integer delay index (scalar)<br>
        `@kmax`: the maximal integer Doppler index (scalar)<br>
        `@kmax_frac`: the maximal float Doppler index (scalar). If set, `kmax` will be forced to set to`floor(kmax_frac)`<br>
        ```matlab
        % matlab
        H_DD = otfs.setChannel("p", 6, "lmax", 11, "kmax", 3);
        ```
        ```python
        # python
        H_DD = otfs.setChannel(p=6, lmax=11, kmax=3);
        ```
    * generate fixed channel<br>
        `@gains:` the channel gains of [(batch_size), p]<br>
        `@delays:` the channel delays of [(batch_size), p]<br>
        `@dopplers`: the channel Doppler shift of [(batch_size), p]<br>
        ```matlab
        % matlab
        H_DD = otfs.setChannel("delays", [0, 1], "Dopplers", [2, 3], "gains", [0.5, 0.5]);
        ```
        ```python
        # python
        # no batch
        H_DD = otfs.setChannel(delays=[0, 1], dopplers=[2, 3], gains=[0.5, 0.5]);
        # using batch
        H_DD = otfs.setChannel(delays=np.tile([0, 1], (batch_size, 1)), dopplers=np.tile([2, 3], (batch_size, 1)), gains=np.tile([0.5, 0.5], (batch_size, 1)));
        ```
    * return value: the channel matrix of [(batch_size), nSubcarNum*nTimeslotNum, nSubcarNum*nTimeslotNum]
* passChannel<br>
    `@No`: a scalar of the linear noise power
* demodulate<br>
    * return value: the vectorized received signal in the delay Doppler domain, a vector of [(batch_size), nSubcarNum*nTimeslotNum]
    ```python
    # matlab & python
    yDD = otfs.demodulate();
    ```
## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require `init` first to load directories.
* TestFractionalDoppler
* TestOTFSAllFunctions: test all functions using random channels and fixed channels
* TestWaveForms
* Viterbo: this code is from `Emanuele Viterbo`. You can download his original code in [this page](https://ecse.monash.edu/staff/eviterbo/OTFS-VTC18/index.html).
