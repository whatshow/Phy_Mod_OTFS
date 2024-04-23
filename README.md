# OTFS Modulation
[![PyPi](https://img.shields.io/badge/PyPi-1.0.3-blue)](https://pypi.org/project/whatshow-phy-mod-otfs/) [![MathWorks](https://img.shields.io/badge/MathWorks-1.0.3-red)](https://mathworks.com/matlabcentral/fileexchange/161136-whatshow_phy_mod_otfs)

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
All codes are uniform in matlab and python in three class.
* OTFSResGrid: this class provides the resource grid
    * OTFSResGrid()<br>
        @in1: 1st input, a scalar for subcarrier number or the content directly<br>
        @in2: only if 1st input is scalar, this input is the `nTimeslotNum`<br>
        ```c, matlab, python
        rg = OTFSResGrid(M, N);     // build a RG (give subcarrier number and timeslot number)
        rg = OTFSResGrid(X_DD);     // build a RG (give X_DD matrix directly)
        ```
    * Set the location of pilots (in the center by default)
        ```c, matlab, python
        rg.pilot2center();  // pilots in the center of the entire OTFS frame
        rg.pilot2zp();      // pilots in the center of zero-padding area
        ```
    * Set the pulse type
        ```c, matlab, python
        rg.setPulse2Ideal();    // use the ideal pulse
        rg.setPulse2Recta();    // use the rectangular pulse
        ```
    * map(): `symbols` is mandatory, `pilots_pow` and `pilots_pow` must use one if you want to use pilots. Other parameters define the length of pilots and the guards area around the pilots. If the guard area is oversize, you are suggested to use full guard on that axis so that the channel estimation result will be more accurate.<br>
        `@symbols`: OTFS symbols<br>
        `@pilots`: a vector of your pilots (if given `pilots_pow` won't be used)<br>
        `@pilots_pow`: pilot power to generate random pilots<br>
        `@pilots_num_delay`: pilots number along the delay(Doppler) axis<br>
        `@pilots_num_doppl`: pilots number along the Doppler(time) axis<br>
        `@guard_delay_full`: full guard on delay (if set true, ignore the number setting)<br>
        `@guard_delay_num_neg`: guard number negatively along the delay(frequency) axis<br>
        `@guard_delay_num_pos`: guard number positively along the delay(frequency) axis<br>
        `@guard_doppl_full`: full guard on Doppler (if set true, ignore the number setting)<br>
        `@guard_doppl_num_neg`: guard number negatively along the Doppler(time) axis<br>
        `@guard_doppl_num_pos`: guard number positively along the Doppler(time) axis
    * demap()<br>
        `@isData`: whether give the data (true by default)<br>
        `@isCE`: whether give the channel estimation result(true by default)<br>
        `@threshold`: the threshold to estimate the channel<br>
        ```c, matlab, python
        [y, his_est, lis_est, kis_est] = rg_rx.demap("threshold", 1e-10);
        y, his_est, lis_est, kis_est = rg_rx.demap("threshold", 1e-10);
        ```
    * clone(): clone this resource grid
    * isZP(): check whether use zero padding or not
    * isPG(): check whether use pilots & guards
    * isInAreaPG()
    
* OTFS: this class provides the entire process of OTFS from Tx to Rx
    * modulate()<br>
    * setChannel(): this methods has 2 kinds of inputs.<br>
    * passChannel()<br>
        `@No`: linear noise power (a scalar) or a given noise vector
        ```c, matlab, python
        otfs.passChannel(No);
        ```
    * demodulate(): return the vectorized received signal in the delay Doppler domain, a vector of [(batch_size), nSubcarNum*nTimeslotNum]
        ```c, matlab, python
        yDD = otfs.demodulate();
        ```
* OTFSDetect: this class provides dedicated OTFS detectors

## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require running `init` in the command window first to load directories.
* `MP_Embed`: examples showing how MP detects symbols using embedded pilots.
* `Tests`
    * `./CE`: examples showing how to estimate channel (using **OTFSResourceGrid**). 
    * `./CE_Detect`: examples showing how to estimate channel (using **OTFSResourceGrid**) and detect symbols. 
    * `test_ce_getcsi`: test `getCSI()` using `sort`
    * `test_ce_detect_*`: test channel estimation and detection together
    * `test_ch_*`: test the channel functions for `OTFS`
    * `test_detect_mp`: test MP detection methods
    * `test_wf` & `test_wf_*`: test waveforms
* TestFractionalDoppler
* `Viterbo_MP_2018`: this code is from `Emanuele Viterbo`. You can download his original code in [this page](https://ecse.monash.edu/staff/eviterbo/OTFS-VTC18/index.html).
* `WaveForm`: this is to observe the waveform of OTFS.