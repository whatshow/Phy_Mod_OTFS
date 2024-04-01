# OTFS Modulation
[![PyPi](https://img.shields.io/badge/PyPi-1.0.1-blue)](https://pypi.org/project/whatshow-phy-mod-otfs/) [![MathWorks](https://img.shields.io/badge/MathWorks-1.0.1-red)](https://mathworks.com/matlabcentral/fileexchange/161136-whatshow_phy_mod_otfs)

This repository is a fundamental toolbox of OTFS modulation crossing `matlab` and `python`. This repositiory is based on papers below:
> Raviteja, P., Phan, K. T., Hong, Y., & Viterbo, E. (2018). Interference cancellation and iterative detection for orthogonal time frequency space modulation. *IEEE transactions on wireless communications, 17(10)*, 6501-6515.

> Raviteja, P., Phan, K. T., & Hong, Y. (2019). Embedded pilot-aided channel estimation for OTFS in delay–Doppler channels. *IEEE transactions on vehicular technology, 68(5)*, 4906-4917.
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
        ```sh, c, matlab, python
        % matlab
        addpath("Modules/OTFS");
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
        `@pilot_loc_type`: pilot location type<br>
        `@batch_size`**(optional)** : the batch size **(only used in python)**.
    
        ```c, matlab, python
        nSubcarNum = 16;
        nTimeslotNum = 7;
        // initialise the class (base)
        otfs = OTFS(nSubcarNum, nTimeslotNum);
        otfs = OTFS(nSubcarNum, nTimeslotNum, batch_size=4); # using the batch (only used in python)
        // initialise the class (base + channel estimation + symbol detection)
        % matlab
        otfs = OTFS(nSubcarNum, nTimeslotNum, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
        # python
        otfs = OTFS(nSubcarNum, nTimeslotNum, pilot_loc_type=OTFS.PILOT_LOC_CENTER);
        ```
    * insertPilotsAndGuards()<br>
        `@pilots`: a vector of your pilots (if given `pilots_pow` won't be used)<br>
        `@pilots_pow`: pilot power to generate random pilots<br>
        `@pilots_num_delay`: pilots number along the delay(Doppler) axis<br>
        `@pilots_num_doppler`: pilots number along the Doppler(time) axis<br>
        `guard_delay_full`: full guard on delay (if set true, ignore the number setting)<br>
        `@guard_delay_num_neg`: guard number negatively along the delay(frequency) axis<br>
        `@guard_delay_num_pos`: guard number positively along the delay(frequency) axis<br>
        `@guard_doppler_full`: full guard on Doppler (if set true, ignore the number setting)<br>
        `@guard_doppler_num_neg`: guard number negatively along the Doppler(time) axis<br>
        `@guard_doppler_num_pos`: guard number positively along the Doppler(time) axis
    
        ```c, matlab, python
        % matlab
        otfs.insertPilotsAndGuards(2, 2, "pilots_pow", 1.2, "guard_delay_num_neg", 1, "guard_delay_num_pos", 1, "guard_doppler_num_neg", 2, "guard_doppler_num_pos", 2);
        # python
        otfs.insertPilotsAndGuards(2, 2, pilots_pow=1.2, guard_delay_num_neg=1, guard_delay_num_pos=1, guard_doppler_num_neg=2, guard_doppler_num_pos=2);
        ```
    * modulate()<br>
        `@symbols`: a vector of `n` symbols as [(batch_size), n], or a matrix in delay Doppler domain [(batch_size), Doppler, delay] or [(batch_size), nTimeslotNum ,nSubcarNum] (the pilots & guards area must be empty)
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
        `@No`: linear noise power (a scalar) or a given noise vector
        ```c, matlab, python
        otfs.passChannel(No);
        ```
    * demodulate(): return the vectorized received signal in the delay Doppler domain, a vector of [(batch_size), nSubcarNum*nTimeslotNum]
        ```c, matlab, python
        yDD = otfs.demodulate();
        ```
    * estimateChannel(): estimate the channel state information based (the pilot type defines the detection method).<br>
        `@threshold:` the threshold of detecting a path
        * The algorithm for SISO single pilot is implemented from<br>
            **III. EMBEDDED CHANNEL ESTIMATION FOR POINT-TO-POINT SISO CASE**<br>
            **C. OTFS With Rectangular Waveforms**
            > Raviteja, P., Phan, K. T., & Hong, Y. (2019). Embedded pilot-aided channel estimation for OTFS in delay–Doppler channels. *IEEE transactions on vehicular technology, 68(5)*, 4906-4917.
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
    * addChannelPath(): add a path to the channel (this does not influence other existing paths)<br>
        `@hi`: the path gain (linear gain)<br>
        `@li`: the delay<br>
        `@ki`: the Doppler shift<br>
    * getChannel(): return the Delay-Doppler domain channel matrix of [(batch_size), n, n]<br>
        `@chan_coef`: the channel gains<br>
        `@delay_taps`: the channel delays<br>
        `@doppler_taps`: the channel Dopplers<br>
        `@only_for_data`: whether the channel is only for data (by default true). If you want to get the entire H_DD when using pilos and/or guards, you should manullay set it to false.
        
        ```c, matlab, python
        H_DD = otfs.getChannel();
        H_DD = otfs.getChannel("chan_coef", chan_coef, "delay_taps", delay_taps, "doppler_taps", doppler_taps, "only_for_data", false); % matlab
        H_DD = otfs.getChannel(chan_coef=chan_coef, delay_taps=delay_taps, doppler_taps=doppler_taps, only_for_data=False); # python
        ```
    * getCSI(): get the channel state information, return [gains, delays, dopplers]<br>
        `@sort_by_gain`: sort axis<br>
        `@sort_by_delay_doppler`: sort axes<br>
        `@sort_by_doppler_delay`: sort axes<br>
        `@sort_ascend`: sort direction<br>
        `@sort_descend`: sort direction
    * getX2DT(): get Tx signal in the Delay Time domain [(batch_size), delay, time]<br>
    * getX2TF(): get Tx signal in the TF domain [(batch_size), nSubcarNum, nTimeslotNum]<br>
    * getX2T(): get Tx signal in the time domain<br>
        `@fft_size`: the size of FFT
    * getYDD(): get Rx received signal in delay Doppler domain<br>

## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require running `init` in the command window first to load directories.
* `Tests`
    * `test_ce_case_01`: test all pilots allocation for **center pilots and reduced guard** (no noise)
    * `test_ce_case_10`: test illegal pilot allocation
    * `test_ce_getcsi`: test `getCSI()` using `sort`
    * `test_ce_detect_*`: test channel estimation and detection together
    * `test_ch_*`: test the channel functions for `OTFS`
    * `test_detect_mp`: test MP detection methods
    * `test_wf` & `test_wf_*`: test waveforms
* TestFractionalDoppler
* `Viterbo_MP_2018`: this code is from `Emanuele Viterbo`. You can download his original code in [this page](https://ecse.monash.edu/staff/eviterbo/OTFS-VTC18/index.html).
* `WaveForm`: this is to observe the waveform of OTFS.