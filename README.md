# OTFS Modulation
[![PyPi](https://img.shields.io/badge/PyPi-2.1.4-blue)](https://pypi.org/project/whatshow-phy-mod-otfs/) [![MathWorks](https://img.shields.io/badge/MathWorks-2.1.4-red)](https://mathworks.com/matlabcentral/fileexchange/161136-whatshow_phy_mod_otfs)

This repository is a fundamental toolbox of OTFS modulation crossing `matlab` and `python`. This repositiory is based on papers below:
> Raviteja, P., Phan, K. T., Hong, Y., & Viterbo, E. (2018). Interference cancellation and iterative detection for orthogonal time frequency space modulation. *IEEE transactions on wireless communications, 17(10)*, 6501-6515.

> Raviteja, P., Phan, K. T., & Hong, Y. (2019). Embedded pilot-aided channel estimation for OTFS in delay–Doppler channels. *IEEE transactions on vehicular technology, 68(5)*, 4906-4917.

## How to install
* Install through `Matlab`
    * `HOME/Add-Ons/Get Add-Ons`: search `whatshow_toolbox` and install it.
    * `HOME/Add-Ons/Get Add-Ons`: search `whatshow_phy_mod_otfs` and install it.
* Install through `pip`
    ```sh
    pip install whatshow-toolbox
    pip install whatshow-phy-mod-otfs
    ```
    * **import this module**
        ```
        from whatshow_phy_mod_otfs import OTFS
        ```

## How to use
All codes are uniform in matlab and python in three class.
![system_model](https://github.com/whatshow/Phy_Mod_OTFS/blob/main/Img/system_model.jpg)
* OTFSResGrid: this class provides the resource grid
    * OTFSResGrid()<br>
        `@in1`: 1st input, a scalar for subcarrier number or the content directly<br>
        `@in2`: only if 1st input is scalar, this input is the `nTimeslotNum`<br>
        `@zp_len(opt)`: zero padding length
        `@batch_size(opt):` the batch size (only used in python)
        ```c, matlab, python
        rg = OTFSResGrid(M, N);     // build a RG (give subcarrier number and timeslot number)
        rg = OTFSResGrid(X_DD);     // build a RG (give X_DD matrix directly)
        // user zero padding
        rg = OTFSResGrid(X_DD, "zp_len", zp_len);
        rg = OTFSResGrid(X_DD, zp_len=zp_len);
        ```
    * Set the pulse type
        ```c, matlab, python
        rg.setPulse2Ideal();    // use the ideal pulse
        rg.setPulse2Recta();    // use the rectangular pulse
        ```
    * Set pilot type
        ```c, matlab, python
        rg.setPilot2Embed();            // use embedded pilots
        rg.setPilot2SuperImposed();     // use superimposed pilots
        ```
    * Set the location of pilots (in the center by default)<br>
        `@pl_len`: pilot length on the delay<br>
        `@pk_len`: pilot length on the doppler<br>
        `@pl1`: pilot location on the delay<br>
        `@pk1`: pilot location on the doppler
        ```c, matlab, python
        rg.pilot2center(pl_len, pk_len);            // pilots in the center of the entire OTFS frame
        rg.pilot2zp(pl_len, pk_len);                // pilots in the center of zero-padding area
        rg.setPilot2Flex(pl_len, pk_len, pl1, pk1); // pilots in user-defined area
        ```
    * setGuard(): set the guard range<br>
        `@in1`: (1,2) negative guard on the delay (3) negative guard on the Doppler<br>
        `@in2`: (1,2) positive guard on the delay (3) positive guard on the Doppler<br>
        `@in3`: (1) negative guard on the Doppler<br>
        `@in4`: (1) positive guard on the Doppler<br>
        `@guard_delay_full(opt)`: full guard on delay (if set true, ignore the number setting)<br>
        `@guard_doppl_full(opt)`: full guard on Doppler (if set true, ignore the number setting)
        ```c, matlab, python
        // set guard range manually
        rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos);
        // full guard
        rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, 'guard_doppl_full', true);
        rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
        ```
    * map(): `symbols` is mandatory, `pilots_pow` and `pilots_pow` must use one if you want to use pilots. Other parameters define the length of pilots and the guards area around the pilots. If the guard area is oversize, you are suggested to use full guard on that axis so that the channel estimation result will be more accurate.<br>
        `@symbols`: OTFS symbols<br>
        `@pilots(opt)`: a vector of your pilots (if given `pilots_pow` won't be used)<br>
        `@pilots_pow(opt)`: pilot power to generate random pilots
    * setAreaCE(): set the channel estimation area manually. **You should only call this function when you disagree with our channel estimation area.**<br>
        `@ce_l_beg`: CE delay beginning<br>
        `@ce_l_end`: CE delay ending<br>
        `@ce_k_beg`: CE Doppler beginning<br>
        `@ce_k_end`: CE Doppler ending<br>
    * demap()<br>
        `@isData(opt)`: whether give the data (true by default)<br>
        `@isCE(opt)`: whether give the channel estimation result(true by default)<br>
        `@threshold(opt)`: the threshold to estimate the channel
        ```c, matlab, python
        [y, his_est, lis_est, kis_est] = rg_rx.demap("threshold", 1e-10);
        y, his_est, lis_est, kis_est = rg_rx.demap(threshold=1e-10);
        ```
    * clone(): clone this resource grid
    * isZP(): return zero padding length
    * isPG(): check whether use pilots & guards
    * Check whether your coordinate is in a certain area<br>
        `@pos_doppl`: the position on the Doppler axis. Not given means the position is for a Doppler-delay vector <br>
        `@pos_delay`: the position on the delay axis for matrix or th position on the Doppler-delay axis for the vector
        ```c, matlab, python
        // is in pilots and guards area
        rg.isInAreaPG(k, l);
        rg.isInAreaPG(n);               // you coordinate is from a vector
        // is in zero padding area
        rg.isInAreaZP(k, l);
        rg.isInAreaZP(n);               // you coordinate is from a vector
        // is in data area
        rg.isInAreaDA(k, l);
        rg.isInAreaDA(n);               // you coordinate is from a vector
        // is in channel estimation area
        rg.isInAreaCE(k, l);
        rg.isInAreaCE(n);               // you coordinate is from a vector
        ```
    * getAreaPG(): get the area of pilots and guards
        ```c, matlab, python
        [pg_num, pg_delay_beg, pg_delay_end, pg_doppl_beg, pg_doppl_end] = rg.getAreaPG();
        pg_num, pg_delay_beg, pg_delay_end, pg_doppl_beg, pg_doppl_end = rg.getAreaPG();
        ```
    * getAreaCE(): get the area of channel estimation
        ```c, matlab, python
        [ce_num, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end] = rg.getAreaCE();
        ce_num, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end = rg.getAreaCE();
        ```
    * check pulse type
        ```c, matlab, python
        rg.isPulseIdeal();  // use ideal pulse or not
        rg.isPulseRecta();  // use rectangula pulse or not
        ```
    * setContent()
        `@content`: a 2D matrix containing pilots, guards and data (if used)
        ```c, matlab, python
        rg.setContent(content);
        ```
    * getContentSize(): get content size
        ```c, matlab, python
        [nSubcarNum, nTimeslotNum] = rg.getContentSize();
        nSubcarNum, nTimeslotNum = rg.getContentSize();
        ```
    * getContent(): get the content<br>
        `@isVector(opt)`: if true, the returned result is a vector
        ```c, matlab, python
        rg.getContent();
        // vectorized results
        rg.getContent("isVector", true);
        rg.getContent(isVector=True);
        ```
    * getContentCE(): get the channel estimation area content (return a matrix)
    * getContentNoCE(): get the content except CE area (return a vector)
    * getContentZeroPG(): get the content of zero PG area (return a vector)
    * getContentZeroCE(): get the content of zero CE area (return a vector)
* OTFS: this class provides the entire process of OTFS from Tx to Rx
    * OTFS()<br>
        `@batch_size(opt)`: the batch size (only used in python)
    * Pulse settings: if you want to input delay-Doppler domain matrix directly into `OTFS`, you need to set the pulse type before `modulate`.
        ```c, matlab, python
        otfs.setPulse2Ideal();  // use ideal pulses
        otfs.setPulse2Recta();  // use rectangular pulses
        ```
    * modulate():modulate (use fast method by default)<br>
        `@in1`: an OTFS resource grid or a 2D matrix [(batch_size), Doppler, delay]<br>
        `@isFast(opt)`: DD domain -> TD domain (no X_TF)
        ```c, matlab, python
        // fast modulate
        otfs.modulate(rg);
        // slow modulate
        otfs.modulate(rg, "isFast", false);
        otfs.modulate(rg, isFast=False);
        ```
    * setChannel(): this methods has multiple kinds of inputs.<br>
        * set a fixed chanel (at least two paths, if you want to add one fixed path, call `setChannelExtra`)
            `@in1->his`: the path gains<br>
            `@in2->lis`: the delays<br>
            `@in3->kis`: the doppler shifts<br>
        * set a random channel (overwritten the channel setting; use Rayleigh fading if not select channel model)
            `@in1->p`: the path number<br>
            `@in2->lmax`: the maxmimal delay index<br>
            `@in3->kmax`: the maximal Doppler index (can be fractional)<br>
            `@force_frac(opt)`: use fractional Doppler (force)<br>
            `@isAWGN(opt)`: use awgn<br>
            `@isRician(opt)`: use Rician fading<br>
            ```c, matlab, python
            otfs.setChannel(p,lmax,kmax);
            // optional inputs
            otfs.setChannel(p,lmax,kmax, "force_frac", true);
            otfs.setChannel(p,lmax,kmax, force_frac=True);
            otfs.setChannel(p,lmax,kmax, "isAWGN", true);
            otfs.setChannel(p,lmax,kmax, isAWGN=True);
            otfs.setChannel(p,lmax,kmax, "isRician", true);
            otfs.setChannel(p,lmax,kmax, isRician=True);
            ```
    * setChannelExtra(): add a path to the channel (this does not influence other existing paths)<br>
        `@hi`: the path gain (linear gain)<br>
        `@li`: the delay<br>
        `@ki`: the Doppler shift<br>
    * passChannel()<br>
        `@No`: linear noise power (a scalar) or a given noise vector
        ```c, matlab, python
        otfs.passChannel(No);
        ```
    * demodulate(): demodulate (use fast method by default)<br>
        `@isFast(opt)`: TD domain -> DD domain (no Y_TF) 
        ```c, matlab, python
        rg = otfs.demodulate(); // if the modulation uses a resource grid
        Y_DD = otfs.demodulate(); // if the modulation uses a delay-Doppler domain matrix
        ```
    * getChannel(): return the channel on delay Doppler domain. If not given `his`, `lis`, `kis`, use the current channel.<br>
        `@his`: the channel gains<br>
        `@lis`: the channel delays<br>
        `@kis`: the channel Dopplers<br>
        `@data_only(opt)`: whether the channel is only for data (by default true). If you want to get the entire H_DD when using pilos and/or guards, you should manullay set it to false.
        ```c, matlab, python
        // the channel only for data
        Hdd = otfs.getChannel("data_only", true);
        Hdd = otfs.getChannel(data_only=True);
        ```
    * getCSI(): get the channel state information
        `@sort_by_gain(opt)`: sort axis, false by defaut<br>
        `@sort_by_delay_doppler(opt)`: sort axes, false by defaut<br>
        `@sort_by_doppler_delay(opt)`: sort axes, false by defaut<br>
        `@descend(opt)`: sort direction, false by defaut
        ```c, matlab, python
        otfs.getCSI("sort_by_delay_doppler", true);
        otfs.getCSI(sort_by_delay_doppler=True);
        ```
    * getXTF(): get the signal in the TF domain
    * getXT(): get the signal in the time domains
        `@fft_size(opt)`: the size of fft
* OTFSDetect: this class provides dedicated OTFS detectors
    * OTFSDetector()<br>
        `@constel`: the constellation (a vector)
        ```c, matlab, python
        od = OTFSDetector(sympool);
        ```
    * Select a detector
        * useMPBase(): set detector types - MP Base<br>
            `@n_ite(opt)`: the iteration number (200 by default)<br>
            `@delta_fra(opt)`: the percentage for taking the values in the current iteration
            ```c, matlab, python
            od.useMPBase();
            od.useMPBase("n_ite", 10, "delta_fra", 0.9);
            od.useMPBase(n_ite=10, delta_fra=0.9);
            ```
    * detect(): return the estimated symbols<br>
        `@y`: the received signal from the resource grid (after demapping) or just a resource grid <br>
        `@csi_info1`: (1) a vector of path gains  (2) a matrix of HDD<br>
        `@csi_info2`: (1) the delay indices       (2) the noise power<br>
        `@csi_info3`: (1) the Doppler indices<br>
        `@csi_info4`: (1) the noise power
        ```c, matlab, python
        xDD_est  = od.detect(rg_rx, his, lis, kis, No);
        xDD_est  = od.detect(rg_rx, Hdd, No);
        ```
    
## Samples
Before running any sample code, please make sure you are at the root path of this repository. Also, Matlab codes require running `init` in the command window first to load directories.
* `Tests`
    * `./Modular_OTFS`: test `OTFS` class
        * `case000`: modular tests on each functions
        * `case001:` test the light mode (integer and fractional Dopplers. ideal and rectangular pulses)    
    * `./Modular_OTFSDetector`: test `OTFSResGrid` class
    * `./Modular_OTFSResGrid`: test `OTFSResGrid` class
    * `./Whole_CE`: examples showing how to estimate channel (using **OTFSResourceGrid**). 
    * `./Whole_CE_Detect`: examples showing how to estimate channel (using **OTFSResourceGrid**) and detect symbols. 
    * `./Whole_Joint`
        * `case001`: ideal/recta pulse + full guards + single pilot; NMSE
        * `case002`: ideal/recta pulse + full guards + single pilot; SER(mp base)
    <!--* `test_ce_getcsi`: test `getCSI()` using `sort`-->
    <!--* `test_ce_detect_*`: test channel estimation and detection together-->
    <!--* `test_ch_*`: test the channel functions for `OTFS`-->
    <!--* `test_detect_mp`: test MP detection methods-->
    <!--* `test_wf` & `test_wf_*`: test waveforms-->
<!--* TestFractionalDoppler-->
* `Viterbo_MP_2018`: this code is from `Emanuele Viterbo`. You can download his original code in [this page](https://ecse.monash.edu/staff/eviterbo/OTFS-VTC18/index.html).
* `Viterbo_MP_Embed`: examples showing how MP detects symbols using embedded pilots.
* `WaveForm`: this is to observe the waveform of OTFS.

## TODO
* `OTFS`
    * getCSI(): python need to sort
    * removeNoDAChannel(): need to remove zero padding area
* `OTFSDetect`
    * `detectMPBase`: support when batch_size is greater than 1