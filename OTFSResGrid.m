% resource grid
classdef OTFSResGrid < handle
    properties(Constant)
        % pulse
        PULSE_NO = 0;
        PULSE_IDEAL = 10;                           % using ideal pulse to estimate the channel
        PULSE_RECTA = 20;                           % using rectangular waveform to estimate the channel
        % pilot types
        PILOT_TYPE_EM = 10;                         % embedded pilots
        PILOT_TYPE_SP = 20;                         % superimposed pilots
        % Pilot locations
        PILOT_LOC_FLEX = 0;                         % flexible location
        PILOT_LOC_CENTER = 10;                      % the pilot is put at the center of frame
        PILOT_LOC_ZP = 20;                          % the pilot is put at the zero padding area
    end
    properties
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        content
        % zero padding
        zp_len = 0;
        % pulse
        pulse_type = OTFSResGrid.PULSE_NO;
        % pilot
        pilot_type = OTFSResGrid.PILOT_TYPE_EM;
        pilots = [];
        p_len = 0;
        pl1 = 0;                                            % 1st (lowest) pilot location in delay axis
        pk1 = 0;                                            % 1st (lowest) pilot location in Doppler axis
        pl_len = 0;                                         % pilots number along the delay(Doppler) axis
        pk_len = 0;                                         % pilots number along the Doppler(time) axis
        pilot_loc_type = OTFSResGrid.PILOT_LOC_CENTER;      % pilot attribute
        % guard
        gl_len_ful = false;
        gl_len_neg = 0;
        gl_len_pos = 0;
        gk_len_ful = false;
        gk_len_neg = 0;
        gk_len_pos = 0;
        % channel estimation area in content
        pg_num = 0;
        pg_delay_beg = NaN;
        pg_delay_end = NaN;
        pg_doppl_beg = NaN;
        pg_doppl_end = NaN;
        % channel estimation area in content
        ce_num = 0;
        ce_delay_beg = NaN;
        ce_delay_end = NaN;
        ce_doppl_beg = NaN;
        ce_doppl_end = NaN;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % general methods
    methods
        %{
        init the resource grid
        @in1:               1st input, a scalar for subcarrier number or the content directly
        @in2:               only if 1st input is scalar, this input is the `nTimeslotNum`
        %}
        function self = OTFSResGrid(in1, varargin)
            % take inputs - nSubcarNum, nTimeslotNum, content & varargin_id_beg
            if isscalar(in1)
                if length(varargin) < 1
                    error("The timeslot number is not given.")
                elseif floor(in1) ~= in1
                    error("The subcarier number can only be an integer.");
                elseif floor(varargin{1}) ~= varargin{1}
                    error("The timeslot number can only be an integer.");
                else
                    self.nSubcarNum = in1;
                    self.nTimeslotNum = varargin{1};
                    self.content = zeros(self.nTimeslotNum, self.nSubcarNum);
                end
            else
                self.content = in1;
                [self.nTimeslotNum, self.nSubcarNum] = size(self.content);
            end
        end
        
        %{
        pulse settings
        %}
        function setPulse2Ideal(self)
            self.pulse_type = self.PULSE_IDEAL;
        end
        function setPulse2Recta(self)
            self.pulse_type = self.PULSE_RECTA;
        end

        %{
        pilot type setting
        %}
        function setPilot2Embed(self)
            self.pilot_type = self.PILOT_TYPE_EM;
        end
        function setPilot2SuperImposed(self)
            self.pilot_type = self.PILOT_TYPE_SP;
        end

        %{
        pilot position setting
        @pl_len:    pilot length on the delay
        @pk_len:    pilot length on the doppler
        @zp_len:    zero padding length
        @pl:        pilot location on the delay
        @pk:        pilot location on the doppler
        %}
        function setPilot2Center(self, pl_len, pk_len)
            self.pilot_loc_type = self.PILOT_LOC_CENTER;
            self.pl_len = pl_len;
            self.pk_len = pk_len;
            self.pl1 = floor((self.nSubcarNum - self.pl_len)/2) + 1;
            self.pk1 = floor((self.nTimeslotNum - self.pk_len)/2) + 1;
            self.validP();
        end
        function setPilot2ZP(self, pl_len, pk_len, zp_len)
            self.pilot_loc_type = self.PILOT_LOC_ZP;
            self.pl_len = pl_len;
            self.pk_len = pk_len;
            self.zp_len = zp_len;
            self.pl1 = floor((self.nSubcarNum - self.pl_len)/2);
            self.pk1 = self.nTimeslotNum - self.zp_len + floor((self.zp_len - self.pk_len)/2);
            self.validP();
        end
        function setPilot2Flex(self, pl_len, pk_len, pl1, pk1)
            self.pilot_loc_type = self.PILOT_LOC_FLEX;
            self.pl_len = pl_len;
            self.pk_len = pk_len;
            self.pl1 = pl1;
            self.pk1 = pk1;
            self.validP();
        end

        %{
        set guards
        @in1:               negative guard on the delay
        @in2:               positive guard on the delay
        @in3:               negative guard on the Doppler
        @in4:               positive guard on the Doppler
        @guard_delay_full:  full guard on delay (if set true, ignore the number setting)
        @guard_doppl_full:  full guard on Doppler (if set true, ignore the number setting)
        %}
        function setGuard(self, varargin)
            % take inputs
            args_len = length(varargin);
            if args_len < 4
                error("The input is not enough.");
            end
            ins = [0,0,0,0];                 % guard lengths on 4 directions
            arg_opt_id_beg = args_len + 1;   % optional arguments beginning id
            for arg_id = 1:args_len
                if isnumeric(varargin{arg_id})
                    ins(arg_id) = varargin{arg_id};
                else
                    arg_opt_id_beg = arg_id;
                    break;
                end
            end

            % take optional inputs
            inPar = inputParser;
            addParameter(inPar, 'guard_delay_full', false, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar, 'guard_doppl_full', false, @(x) isscalar(x)&&islogical(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{arg_opt_id_beg:end});
            self.gl_len_ful = inPar.Results.guard_delay_full;
            self.gk_len_ful = inPar.Results.guard_doppl_full;
            % input check - full guard
            if self.gl_len_ful
                ins(1) = floor((self.nSubcarNum - self.pl_len)/2);
                ins(2) = self.nSubcarNum - self.pl_len - ins(1);               
            end
            if self.gk_len_ful
                ins(3) = floor((self.nTimeslotNum - self.pk_len)/2);
                ins(4) = self.nTimeslotNum - self.pk_len - ins(3);
            end
            % input check - guard - integers only
            if floor(ins(1)) ~= ins(1) || floor(ins(2)) ~= ins(2)
                error("Guard number along the delay axis must be integers.");
            end
            if floor(ins(3)) ~= ins(3) || floor(ins(4)) ~= ins(4)
                error("Guard number along the Doppler axis must be integers.");
            end
            % input check - guard - no negative
            if ins(1) < 0 || ins(2) < 0
                error("Guard number along the delay axis must be non-negative.");
            end
            if ins(3) < 0 || ins(4) < 0
                error("Guard number along the Doppler axis must be non-negative.");
            end           
            % take inputs
            self.gl_len_neg = ins(1);
            self.gl_len_pos = ins(2);
            self.gk_len_neg = ins(3);
            self.gk_len_pos = ins(4);
        end

        %{
        map
        @symbols:                   OTFS symbols
        @pilots:                    a vector of your pilots (if given `pilots_pow` won't be used)
        @pilots_pow:                pilot power to generate random pilots
        %}
        function map(self, symbols, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar, 'pilots', self.pilots, @(x) isvector(x)&&isnumeric(x));
            addParameter(inPar, 'pilots_pow', NaN, @(x) isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{:});
            % optional inputs - assign
            self.pilots         = inPar.Results.pilots;
            pilots_pow          = inPar.Results.pilots_pow;
            
            self.validP();
            self.calcAreaPGCE();
            % insert
            self.insertDA(symbols);
            self.insertP(pilots_pow);
        end

        %{
        set the channel estimate area (CE area is set in `map`, if you want to use your own area, call this)
        @ce_l_beg: CE delay beginning
        @ce_l_end: CE delay ending
        @ce_k_beg: CE Doppler beginning
        @ce_k_end: CE Doppler ending
        %}
        function setAreaCE(self, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end)
            if ce_delay_beg < 1 || ce_delay_beg > self.nSubcarNum
                error("CE delay beginning overflows.");
            elseif ce_delay_end < 1 || ce_delay_end > self.nSubcarNum
                error("CE delay ending overflows.");
            elseif ce_delay_beg > ce_delay_end
                error("CE delay beginning is after CE delay ending.");
            end
            if ce_doppl_beg < 1 || ce_doppl_beg > self.nTimeslotNum
                error("CE Doppler beginning overflows.");
            elseif ce_doppl_end < 1 || ce_doppl_end > self.nTimeslotNum
                error("CE Doppler ending overflows.");
            elseif ce_doppl_beg > ce_doppl_end
                error("CE Doppler beginning is after CE delay ending.");
            end
            self.ce_delay_beg = ce_delay_beg;
            self.ce_delay_end = ce_delay_end;
            self.ce_doppl_beg = ce_doppl_beg;
            self.ce_doppl_end = ce_doppl_end;
            self.ce_num = (ce_delay_end - ce_delay_beg + 1)*(ce_doppl_end - ce_doppl_beg + 1);
        end

        %{
        demap
        @isData:        whether give the data
        @isCE:          whether give the channel estimation result
        @threshold:     the threshold to estimate the channel
        %}
        function [y, his, lis, kis] = demap(self, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar, 'isData', true, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar, 'isCE', true, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar, 'threshold', 0, @(x) isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{:});
            % optional inputs - assign
            isData = inPar.Results.isData;
            isCE = inPar.Results.isCE;
            threshold = inPar.Results.threshold;
            % input check          
            % input check - threshold
            if threshold < 0
                error("The threshould must be non-negative.");
            end
            % input check - pulse_type
            if self.pulse_type == self.PULSE_NO
                error("The pulse type has to be set before demapping.");
            end

            % y
            y = NaN;
            if isData
                y = self.getContentNoCE();
            end
            % Hest
            his = NaN;
            lis = NaN;
            kis = NaN;
            % need to estimate the channel & have pilots
            if isCE && self.isPG()
                if self.p_len == 1
                    [his, lis, kis] = self.estimateChannel4SingPilot(threshold);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % issers, getters, setters
    methods
        %{
        clone
        %}
        function rg = clone(self)
            rg = OTFSResGrid(self.content);
            rg.zp_len = self.zp_len;
            rg.pulse_type = self.pulse_type;
            rg.pilots = self.pilots;
            rg.p_len = self.p_len;
            rg.pl_len = self.pl_len;
            rg.pk_len = self.pk_len;
            rg.pl1 = self.pl1;
            rg.pk1 = self.pk1;
            rg.pilot_loc_type = self.pilot_loc_type;
            rg.pg_num = self.pg_num;
            rg.pg_delay_beg = self.pg_delay_beg;
            rg.pg_delay_end = self.pg_delay_end;
            rg.pg_doppl_beg = self.pg_doppl_beg;
            rg.pg_doppl_end = self.pg_doppl_end;
            rg.ce_num = self.ce_num;
            rg.ce_delay_beg = self.ce_delay_beg;
            rg.ce_delay_end = self.ce_delay_end;
            rg.ce_doppl_beg = self.ce_doppl_beg;
            rg.ce_doppl_end = self.ce_doppl_end;
        end

        %{
        return zero padding length
        %}
        function zp_len = isZP(self)
            zp_len = self.zp_len;
        end

        %{
        check whether use pilots & guards
        %}
        function is_pg = isPG(self)
            % check whether pilots is assigned or not
            is_pg = ~isempty(self.pilots);
            % check whether pilot and guard area is calculated
            is_pg = is_pg && self.pg_num > 0;
            % check whether CE area is calulated
            is_pg = is_pg && self.ce_num > 0;
        end

        %{
        check whether the current position is in channel estimation area
        @pos_doppl:     the position on the Doppler axis. Not given means the position is for a Doppler-delay vector
        @pos_delay:     the position on the delay axis for matrix or th position on the Doppler-delay axis for the vector
        %}
        function is_in = isInAreaPG(self, pos_doppl, varargin)
            is_in = self.isInArea(10, pos_doppl, varargin{:});
        end
        function is_in = isInAreaZP(self, pos_doppl, varargin)
            is_in = self.isInArea(11, pos_doppl, varargin{:});
        end
        function is_in = isInAreaDA(self, pos_doppl, varargin)
            is_in = self.isInArea(12, pos_doppl, varargin{:});
        end
        function is_in = isInAreaCE(self, pos_doppl, varargin)
            is_in = self.isInArea(20, pos_doppl, varargin{:});
        end

        %{
        get the area of pilots and guards
        %}
        function [pg_num, pg_delay_beg, pg_delay_end, pg_doppl_beg, pg_doppl_end] = getAreaPG(self)
            pg_num = self.pg_num;
            pg_delay_beg = self.pg_delay_beg;
            pg_delay_end = self.pg_delay_end;
            pg_doppl_beg = self.pg_doppl_beg;
            pg_doppl_end = self.pg_doppl_end;
        end

        %{
        get the area of channel estimation
        %}
        function [ce_num, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end] = getAreaCE(self)
            ce_num = self.ce_num;
            ce_delay_beg = self.ce_delay_beg;
            ce_delay_end = self.ce_delay_end;
            ce_doppl_beg = self.ce_doppl_beg;
            ce_doppl_end = self.ce_doppl_end;
        end

        %{
        check pulse type
        %}
        function is_pulse_type = isPulseIdeal(self)
            is_pulse_type = self.pulse_type == self.PULSE_IDEAL;
        end
        function is_pulse_type = isPulseRecta(self)
            is_pulse_type = self.pulse_type == self.PULSE_RECTA;
        end
        
        %{
        set content
        %}
        function setContent(self, content)
            self.content = content;
        end

        %{
        get content size
        %}
        function [nSubcarNum, nTimeslotNum] = getContentSize(self)
            nSubcarNum = self.nSubcarNum;
            nTimeslotNum = self.nTimeslotNum;
        end

        %{
        get content
        @isVector: if true, the returned result is a vector
        %}
        function content = getContent(self, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar, 'isVector', false, @(x) isscalar(x)&&islogical(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{:});
            % take inputs
            isvector = inPar.Results.isVector;
            % return
            if isvector
                content = self.content.';
                content = content(:);
            else
                content = self.content;
            end
        end

        %{
        get content - CE
        %}
        function ce_area = getContentCE(self)
            if self.ce_num == 0
                ce_area = [];
            else
                ce_area = self.content(self.ce_doppl_beg:self.ce_doppl_end, self.ce_delay_beg:self.ce_delay_end);
            end
        end
        
        %{
        get content - except CE area (return a vector)
        %}
        function data = getContentNoCE(self)
            % get data
            if self.ce_num == 0
                data = self.content.';
                data = data(:);
            else
                data_num = self.nSubcarNum*self.nTimeslotNum - self.ce_num;
                data = zeros(data_num, 1);
                data_id = 1;
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if ~self.isInAreaCE(doppl_id, delay_id)
                            data(data_id) = self.content(doppl_id, delay_id);
                            data_id = data_id + 1;
                        end
                    end
                end
            end
        end

        %{
        get content - zero PG area (return a vector)
        %}
        function data = getContentZeroPG(self, varargin)
            % get data
            data = self.content;
            % zero the CE area
            if self.ce_num > 0
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if self.isInAreaPG(doppl_id, delay_id)
                            data(doppl_id, delay_id) = 0;
                        end
                    end
                end
            end
        end

        %{
        get content - zero CE area (return a vector)
        %}
        function data = getContentZeroCE(self, varargin)
            % get data
            data = self.content;
            % zero the CE area
            if self.ce_num > 0
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if self.isInAreaCE(doppl_id, delay_id)
                            data(doppl_id, delay_id) = 0;
                        end
                    end
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % private methods
    methods(Access=private)
        %{
        validate pilot settings
        %}
        function validP(self)
            % 1st pilot coordinates & % pilot length
            if self.pk_len > 0 && self.pl_len > 0
                if self.pk1 < 1 || self.pk1 > self.nTimeslotNum
                    error("Pilot 1st location on the delay axis overflows.");
                end
                if self.pl1 < 1 || self.pl1 > self.nSubcarNum
                    error("Pilot 1st location on the Doppler axis overflows.");
                end
                if self.pk1 + self.pk_len - 1 > self.nTimeslotNum
                    error("Pilot length on the delay axis overflows.");
                end
                if self.pl1 + self.pl_len - 1 > self.nSubcarNum
                    error("Pilot length on the Doppler axis overflows.");
                end
            end
            % zero padding
            if self.zp_len < 0 || self.zp_len > self.nSubcarNum
                error("Zero Padding length cannot be negative or over subcarrier number.");
            end
        end

        %{
        calculate PG & CE area
        %}
        function calcAreaPGCE(self)
            % overflow check
            if self.pl1 - self.gl_len_neg <= 0
                error("The guard (neg) on delay axis overflows.");
            end
            if (self.pl1+self.pl_len-1) + self.gl_len_pos > self.nSubcarNum
                error("The guard (pos) on delay axis overflows.");
            end
            if self.pk1 - self.gk_len_neg <= 0
                error("The guard (neg) on Doppler axis overflows.");
            end
            if (self.pk1+self.pk_len-1) + self.gk_len_pos > self.nTimeslotNum
                error("The guard (pos) on Doppler axis overflows.");
            end
            % calculate area
            if self.pl_len > 0 && self.pk_len > 0
                % calculate PG area
                if self.pilot_type == self.PILOT_TYPE_EM
                    % PG area only exist when using embedded pilots
                    self.pg_num = (self.pl_len+self.gl_len_neg+self.gl_len_pos)*(self.pk_len+self.gk_len_neg+self.gk_len_pos);
                    self.pg_delay_beg = self.pl1 - self.gl_len_neg;
                    self.pg_delay_end = self.pl1 + self.pl_len - 1 + self.gl_len_pos;
                    self.pg_doppl_beg = self.pk1 - self.gk_len_neg;
                    self.pg_doppl_end = self.pk1 + self.pk_len - 1 + self.gk_len_pos;
                end
                
                % calulate channel estimate area
                if self.gl_len_ful
                    self.ce_delay_beg = 1;
                    self.ce_delay_end = self.nSubcarNum;
                else
                    self.ce_delay_beg = self.pg_delay_beg + self.gl_len_neg;
                    self.ce_delay_end = self.pg_delay_end;
                end
                if self.gk_len_ful
                    self.ce_doppl_beg = 1;
                    self.ce_doppl_end = self.nTimeslotNum;
                else
                    self.ce_doppl_beg = self.pg_doppl_beg + floor(self.gk_len_neg/2);
                    self.ce_doppl_end = self.pg_doppl_end - floor(self.gk_len_pos/2);
                end
                self.ce_num = (self.ce_delay_end - self.ce_delay_beg + 1)*(self.ce_doppl_end - self.ce_doppl_beg + 1);
            end
        end
        
        %{
        insert data
        @symbols: symbols to map a vector
        %}
        function insertDA(self, symbols)
            % input check
            if self.pilot_type == self.PILOT_TYPE_SP
                data_num = self.nTimeslotNum*self.nSubcarNum - self.zp_len*self.nTimeslotNum;
            else
                data_num = self.nTimeslotNum*self.nSubcarNum - self.pg_num - self.zp_len*self.nTimeslotNum;
            end
            % input check - symbols
            if ~isvector(symbols)
                error("The transmission symbol must be a vector");
            elseif length(symbols) ~= data_num
                error("The transmission symbol number must be %d", data_num);
            end
            
            % modulate
            % reshape(rowwise) to [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
            %TODO: fill data when using zero padding
            symbols = symbols(:);
            if data_num == self.nSubcarNum*self.nTimeslotNum
                self.content = transpose(reshape(symbols, self.nSubcarNum, self.nTimeslotNum));
            elseif data_num == self.nSubcarNum*self.nTimeslotNum - self.zp_len*self.nTimeslotNum
                self.content = [transpose(reshape(symbols, (self.nSubcarNum - self.zp_len), self.nTimeslotNum)); zeros(self.zp_len, self.nTimeslotNum)];
            else
                symbols_id = 1;
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if self.isInAreaDA(doppl_id, delay_id)
                            self.content(doppl_id, delay_id) = symbols(symbols_id);
                            symbols_id = symbols_id + 1;
                        end
                    end
                end
                assert(symbols_id - 1 == data_num);
            end
        end

        %{
        insert pilots
        @pilots_pow: pilot power
        %}
        function insertP(self, pilots_pow)
            % input check
            self.p_len = length(self.pilots);
            % input check - pilots
            if ~isempty(self.pilots)
                if self.p_len ~= self.pl_len*self.pk_len
                    error("The manual pilot input do not have the required number.");
                elseif self.p_len > self.nTimeslotNum*self.nSubcarNum
                    error("The manual pilot input overflows (over the OTFS frame size).");
                end
            end
            % input check - pilot power
            if self.p_len == 0 && isnan(pilots_pow)
                error("The pilots (linear) power is required while no manual pilot input.");
            end
            % initiate pilots if empty
            if isempty(self.pilots)
                self.p_len = self.pl_len*self.pk_len;
                self.pilots = sqrt(pilots_pow/2)*(1+1j)*ones(self.p_len, 1);
            end
            % allocate pilots
            if self.p_len ~= 0
                % allocate pilots
                self.content(self.pk1:self.pk1+self.pk_len-1, self.pl1:self.pl1+self.pl_len-1) = self.content(self.pk1:self.pk1+self.pk_len-1, self.pl1:self.pl1+self.pl_len-1) + transpose(reshape(self.pilots, self.pl_len, self.pk_len));
            end
        end

        %{
        decide whether the current location is in CE area
        @tag:           10->PG, 11->ZP, 12->Data, 20->CE
        @pos_doppl:     the position on the Doppler axis.
        @pos_delay:     the delay position. Not given means the position is for a Doppler-delay vector
        %}
        function is_in = isInArea(self, tag, pos_doppl, varargin)
            % input check
             if isempty(varargin)
                if pos_doppl <= 0 || pos_doppl > self.nSubcarNum*self.nTimeslotNum
                    error("The vector position is out of the OTFS size.");
                end
            else
                if pos_doppl <= 0 || pos_doppl > self.nTimeslotNum
                    error("The Doppler position is out of the timeslot number.");
                end
                if varargin{1} <= 0 || varargin{1} > self.nSubcarNum
                    error("The delay position is out of the subcarrier number.");
                end
            end
            
            % recalculate the position
            if isempty(varargin)
                pos_id = pos_doppl;
                pos_doppl = pos_doppl/self.nSubcarNum;
                if pos_doppl == floor(pos_doppl)
                    pos_doppl = floor(pos_doppl);
                else
                    pos_doppl = floor(pos_doppl) + 1;
                end
                delay_pos = pos_id - (pos_doppl-1)*self.nSubcarNum;
            else
                delay_pos = varargin{1};
            end
            % decide
            if tag == 10
                is_in = delay_pos >= self.pg_delay_beg && delay_pos <= self.pg_delay_end && pos_doppl >= self.pg_doppl_beg && pos_doppl <= self.pg_doppl_end;
            elseif tag == 11
                is_in = self.zp_len > 0 && delay_pos > self.nSubcarNum - self.zp_len;
            elseif tag == 12
                if self.pg_num == 0
                    is_in = true;
                else
                    is_in = delay_pos < self.pg_delay_beg || delay_pos > self.pg_delay_end || pos_doppl < self.pg_doppl_beg || pos_doppl > self.pg_doppl_end;
                end
                if self.zp_len > 0 
                    is_in = is_in && delay_pos <= self.nSubcarNum - self.zp_len;
                end
            elseif tag == 20
                is_in = delay_pos >= self.ce_delay_beg && delay_pos <= self.ce_delay_end && pos_doppl >= self.ce_doppl_beg && pos_doppl <= self.ce_doppl_end;
            end
        end

        %{
        % estimated the channel using a sngle pilot
        % @threshold: the threshold to detect a path
        %}
        function [his, lis, kis] = estimateChannel4SingPilot(self, threshold)
            % input check
            if self.p_len ~= 1
                error("There should be only 1 pilot.");
            end
            
            % estimate the channel
            his = [];
            lis = [];
            kis = [];
            for delay_id = self.ce_delay_beg:self.ce_delay_end
                for doppl_id = self.ce_doppl_beg:self.ce_doppl_end
                    pss_y = self.content(doppl_id, delay_id);
                    if abs(pss_y) > threshold
                        li = delay_id - self.pl1;
                        ki = doppl_id - self.pk1;
                        if self.pulse_type == self.PULSE_IDEAL
                            pss_beta = exp(-2j*pi*li*ki/self.nSubcarNum/self.nTimeslotNum);
                        elseif self.pulse_type == self.PULSE_RECTA
                            pss_beta = exp(2j*pi*(self.pl1 - 1)*ki/self.nSubcarNum/self.nTimeslotNum);
                        end
                        hi = pss_y/self.pilots/pss_beta;
                        his(end+1) = hi;
                        lis(end+1) = li;
                        kis(end+1) = ki;
                    end
                end
            end
        end
    end
end