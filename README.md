# OTFS modulation
This repositiory is a foundamental toolbox of OTFS modulation crossing `matlab` and `python`.

## How to install
Currently, we offer three options to install this tool.
* Install through `Matlab Add-Ons`

    Go to [releases](https://github.com/whatshow/Phy_Mod_OTFS/releases) to download the file `whatshow_phy_mod_otfs.mltbx` in the latest release to install.<br>
    The class `OTFS` will be automatically added in Matlab. (Don't forget to restart Matlab after installation).
* Install through `pip`
    ```sh
    pip install whatshow-phy-mod-otfs
    ```
    * **import this module**
        ```
        from whatshow_phy_mod_otfs import OTFS
        ```
* Install through git under another local repositiory
    ```sh
    git submodule add git@github.com:whatshow/Phy_Mod_OTFS.git Modules/Mod_OTFS
    ```
    Now, you can see a folder `Modules` with `Mod_OTFS` inside
    * **Import this model**
        ```
        % matlab
        addpath("Modules/OTFS");
        ```
        ```python
        # python
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
        `@kmax`: the maximal integer Doppler index (scalar). If set to a float number, we use fractional Doppler<br>
        ```matlab
        % matlab
        otfs.setChannel("p", 6, "lmax", 11, "kmax", 3);
        ```
        ```python
        # python
        otfs.setChannel(p=6, lmax=11, kmax=3);
        ```
    * generate fixed channel<br>
        `@gains:` the channel gains of [(batch_size), p]<br>
        `@delays:` the channel delays of [(batch_size), p]<br>
        `@dopplers`: the channel Doppler shift of [(batch_size), p]<br>
        ```matlab
        % matlab
        otfs.setChannel("delays", [0, 1], "Dopplers", [2, 3], "gains", [0.5, 0.5]);
        ```
        ```python
        # python
        # no batch
        otfs.setChannel(delays=[0, 1], dopplers=[2, 3], gains=[0.5, 0.5]);
        # using batch
        otfs.setChannel(delays=np.tile([0, 1], (batch_size, 1)), dopplers=np.tile([2, 3], (batch_size, 1)), gains=np.tile([0.5, 0.5], (batch_size, 1)));
        ```
* getChannel<br>
    return the Delay-Doppler domain channel matrix of [(batch_size), nSubcarNum*nTimeslotNum, nSubcarNum*nTimeslotNum]
    ```python
    H_DD = otfs.getChannel();
    ```
* passChannel<br>
    `@No`: a scalar of the linear noise power
    ```python
    # matlab & python
    otfs.passChannel(No);
    ```
* demodulate<br>
    * return value: the vectorized received signal in the delay Doppler domain, a vector of [(batch_size), nSubcarNum*nTimeslotNum]
    ```python
    # matlab & python
    yDD = otfs.demodulate();
    ```
## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require running `init` in the command window first to load directories.
* `TestDetect`: test all kinds of OTFS detectors
* TestFractionalDoppler
* TestOTFSAllFunctions: test all functions using random channels and fixed channels
* TestWaveForms
* `Viterbo_MP_2018`: this code is from `Emanuele Viterbo`. You can download his original code in [this page](https://ecse.monash.edu/staff/eviterbo/OTFS-VTC18/index.html).
