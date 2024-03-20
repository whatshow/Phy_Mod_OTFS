# OTFS modulation
[![PyPi](https://img.shields.io/badge/PyPi-1.0.1-blue)](https://pypi.org/project/whatshow-phy-mod-otfs/) [![MathWorks](https://img.shields.io/badge/MathWorks-1.0.1-red)](https://mathworks.com/matlabcentral/fileexchange/161136-whatshow_phy_mod_otfs)

This repository is a fundamental toolbox of OTFS modulation crossing `matlab` and `python`.

## How to install
Currently, we offer three options to install this tool.
* Install through `Matlab Add-Ons`
    * Install through Matlab `Get Add-Ons`: search `whatshow_phy_mod_otfs` and install it.
    * Install through `.mltbx`: Go to [releases](https://github.com/whatshow/Phy_Mod_OTFS/releases) to download the file `whatshow_phy_mod_otfs.mltbx` in the latest release to install.
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
All codes are uniform in matlab and python as a class.
* This section illustrate the methods of the class following the process from Tx to Rx.
    * OTFS()<br>
        `@nSubcarNum`: the subcarrier number<br>
        `@nTimeslotNum`: the timeslot number<br>
        `pilot_type`: pilot type<br>
        `@pilot_loc_type`: pilot location type<br>
        `@guard_type`: guard type<br>
        `@batch_size`**(optional)** : the batch size **(only used in python)**.
    
        ```c, matlab, python
        nSubcarNum = 16;
        nTimeslotNum = 7;
        // initialise the class (base)
        otfs = OTFS(nSubcarNum, nTimeslotNum);
        otfs = OTFS(nSubcarNum, nTimeslotNum, batch_size=4); # using the batch (only used in python)
        // initialise the class (base + channel estimation + symbol detection)
        % matlab
        otfs = OTFS(nSubcarNum, nTimeslotNum, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
        # python
        otfs = OTFS(nSubcarNum, nTimeslotNum, pilot_type=OTFS.PILOT_SINGLE_SISO, pilot_loc_type=OTFS.PILOT_LOC_CENTER, guard_type=OTFS.GUARD_REDUCED);
        ```
    * insertPilotsAndGuards()<br>
        `@pilots`: a vector of your pilots (if given `pilots_pow` won't be used)<br>
        `@pilots_pow`: pilot power to generate random pilots<br>
        `@pilots_num_delay`: pilots number along the delay(Doppler) axis<br>
        `@pilots_num_doppler`: pilots number along the Doppler(time) axis<br>
        `@guard_delay_num_neg`: guard number negatively along the delay(frequency) axis<br>
        `@guard_delay_num_pos`: guard number positively along the delay(frequency) axis<br>
        `@guard_doppler_num_neg`: guard number negatively along the Doppler(time) axis<br>
        `@guard_doppler_num_pos`: guard number positively along the Doppler(time) axis
    
        ```c, matlab, python
        SNR_p = 30; % dB
        pil_pow = 10^(SNR_p/10);
        lmax = 1;
        kmax = 1;
        pilots_num_delay = 2;
        pilots_num_doppler = 2;
        guard_delay_num_neg = lmax;
        guard_delay_num_pos = lmax;
        guard_doppl_num_neg = kmax*2;
        guard_doppl_num_pos = kmax*2;
        otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
        ```
    * modulate()<br>
        `@symbols`: a vector of [(batch_size), nSubcarNum*nTimeslotNum], or a matrix in delay Doppler domain [(batch_size), Doppler, delay] or [(batch_size), nTimeslotNum ,nSubcarNum]
        ```c, matlab, python
        otfs.modulate(x_origin);
        ```
    * setChannel(): this methods has 2 kinds of inputs.<br>
        * generate random channel<br>
            If we use `batch`, the channel will be different for each sample in the batch.<br>
            `@p`: the path number (scalar)<br>
            `@lmax`: the maximal integer delay index (scalar)<br>
            `@kmax`: the maximal integer Doppler index (scalar). If set to a float number, we use fractional Doppler
            
            ```c, matlab, python
            % matlab
            otfs.setChannel("p", 6, "lmax", 11, "kmax", 3);
            # python
            otfs.setChannel(p=6, lmax=11, kmax=3);
            ```
        * generate fixed channel<br>
            `@gains:` the channel gains of [(batch_size), p]<br>
            `@delays:` the channel delays of [(batch_size), p]<br>
            `@dopplers`: the channel Doppler shift of [(batch_size), p]
            
            ```c, matlab, python
            % matlab
            otfs.setChannel("delays", [0, 1], "Dopplers", [2, 3], "gains", [0.5, 0.5]);
            # python
            otfs.setChannel(delays=[0, 1], dopplers=[2, 3], gains=[0.5, 0.5]);
            otfs.setChannel(delays=np.tile([0, 1], (batch_size, 1)), dopplers=np.tile([2, 3], (batch_size, 1)), gains=np.tile([0.5, 0.5], (batch_size, 1))); # using batch
            ```
    * passChannel()<br>
        `@No`: a scalar of the linear noise power
        ```c, matlab, python
        otfs.passChannel(No);
        ```
    * demodulate(): return the vectorized received signal in the delay Doppler domain, a vector of [(batch_size), nSubcarNum*nTimeslotNum]
        ```c, matlab, python
        yDD = otfs.demodulate();
        ```
    * detect(): detect the symbols from received signal<br>
        `@detect_type`: detect type<br>
        `@csi_type`: detect CSI type<br>
        `@No`: estimated noise power (linear)<br>
        `@constellation`: the constellation (a vector)<br>
        `@chan_coef`: the estimated channel gains<br>
        `@delay_taps`: the estimated channel delays<br>
        `@doppler_taps`: the estimated channel<br>
        
        ```c, matlab, python
        // perfect CSI
        otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_PERFECT, No, sympool);
        // estimated CSI
        otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_CE, No, sympool);
        // input CSI
        % matlab
        otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_IN, No, sympool, "chan_coef", chan_coef, "delay_taps", delay_taps, "doppler_taps", doppler_taps);
        # python
        otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_IN, No, sympool, chan_coef=chan_coef, delay_taps=delay_taps, doppler_taps=doppler_taps);
        ```
* This section illustrate support methods to expose attributes
    * getChannel(): return the Delay-Doppler domain channel matrix of [(batch_size), nSubcarNum*nTimeslotNum, nSubcarNum*nTimeslotNum]
        ```python
        H_DD = otfs.getChannel();
        ```
    * addChannelPath(): add a path to the channel (this does not influence other existing paths)<br>
        `@hi`: the path gain (linear gain)<br>
        `@li`: the delay<br>
        `@ki`: the Doppler shift<br>
    * getChannelGains(): get the channel gains<br>
    * getChannelDelays(): get the channel delays<br>
    * getChannelDopplers(): get the channel dopplers<br>
    * getX2DT(): get Tx signal in the Delay Time domain [(batch_size), delay, time]<br>
    * getX2TF(): get Tx signal in the TF domain [(batch_size), nSubcarNum, nTimeslotNum]<br>
    * getX2T(): get Tx signal in the time domain<br>
        `@fft_size`: the size of FFT
    * getYDD(): get Rx received signal in delay Doppler domain<br>

## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require running `init` in the command window first to load directories.
* `Tests`
    * `test_ce_case_01`: test all pilots allocation for **center pilots and reduced guard** (no noise)
    * `test_detect_mp`: test MP detection methods
* TestFractionalDoppler
* TestOTFSAllFunctions: test all functions using random channels and fixed channels
* TestWaveForms
* `Viterbo_MP_2018`: this code is from `Emanuele Viterbo`. You can download his original code in [this page](https://ecse.monash.edu/staff/eviterbo/OTFS-VTC18/index.html).
