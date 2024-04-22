% resource grid
classdef OTFSResGrid < handle
    properties(Constant)
        % pulse
        PULSE_NO = 0;
        PULSE_IDEAL = 10;                           % using ideal pulse to estimate the channel
        PULSE_RECTA = 20;                           % using rectangular waveform to estimate the channel
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
        pilots = [];
        pilots_len = 0;
        pilots_num_delay = 0;                       % pilots number along the delay(Doppler) axis
        pilots_num_doppl = 0;                       % pilots number along the Doppler(time) axis
        pilot_loc_delay_1st = 0;                    % 1st (lowest) pilot location in delay axis
        pilot_loc_doppl_1st = 0;                    % 1st (lowest) pilot location in Doppler axis
        % pilot attribute
        pilot_loc_type = OTFSResGrid.PILOT_LOC_CENTER;
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
        pilot position setting
        @zp_len: zero padding length
        @pl:    pilot location on the delay
        @pk:    pilot location on the doppler
        %}
        function pilot2Center(self)
            self.pilot_loc_type = self.PILOT_LOC_CENTER;
        end
        function pilot2ZP(self, zp_len)
            self.pilot_loc_type = self.PILOT_LOC_ZP;
            self.zp_len = zp_len;
            if self.zp_len < 0 || self.zp_len > self.nSubcarNum
                error("Zero Padding length cannot be negative or over subcarrier number.");
            end
        end
        function pilot2Flex(self, pl, pk)
            self.pilot_loc_type = self.PILOT_LOC_FLEX;
            self.pilot_loc_delay_1st = pl;
            self.pilot_loc_doppl_1st = pk;
        end

        %{
        map
        @symbols:                   OTFS symbols
        @pilots:                    a vector of your pilots (if given `pilots_pow` won't be used)
        @pilots_pow:                pilot power to generate random pilots
        @pilots_num_delay:          pilots number along the delay(Doppler) axis
        @pilots_num_doppl:          pilots number along the Doppler(time) axis
        @guard_delay_full:          full guard on delay (if set true, ignore the number setting)
        @guard_delay_num_neg:       guard number negatively along the delay(frequency) axis
        @guard_delay_num_pos:       guard number positively along the delay(frequency) axis
        @guard_doppl_full:          full guard on Doppler (if set true, ignore the number setting)
        @guard_doppl_num_neg:       guard number negatively along the Doppler(time) axis
        @guard_doppl_num_pos:       guard number positively along the Doppler(time) axis
        %}
        function map(self, symbols, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar, 'pilots', self.pilots, @(x) isvector(x)&&isnumeric(x));
            addParameter(inPar, 'pilots_pow', NaN, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'pilots_num_delay', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'pilots_num_doppl', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_delay_full', false, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar, 'guard_doppl_full', false, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar, 'guard_delay_num_neg', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_delay_num_pos', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_doppl_num_neg', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_doppl_num_pos', 0, @(x) isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{:});
            % optional inputs - assign
            self.pilots              = inPar.Results.pilots;
            self.pilots_num_delay    = inPar.Results.pilots_num_delay;
            self.pilots_num_doppl    = inPar.Results.pilots_num_doppl;
            pilots_pow          = inPar.Results.pilots_pow;
            guard_delay_full    = inPar.Results.guard_delay_full;
            guard_doppl_full    = inPar.Results.guard_doppl_full;
            guard_delay_num_neg = inPar.Results.guard_delay_num_neg;
            guard_delay_num_pos = inPar.Results.guard_delay_num_pos;
            guard_doppl_num_neg = inPar.Results.guard_doppl_num_neg;
            guard_doppl_num_pos = inPar.Results.guard_doppl_num_pos;
            
            % insert pilots and guards
            self.insertP(pilots_pow);
            self.insertG(guard_delay_full, guard_doppl_full, guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos);
            % insert data
            self.insertDA(symbols);
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
                if self.pilots_len == 1
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
            rg.pilots_len = self.pilots_len;
            rg.pilots_num_delay = self.pilots_num_delay;
            rg.pilots_num_doppl = self.pilots_num_doppl;
            rg.pilot_loc_delay_1st = self.pilot_loc_delay_1st;
            rg.pilot_loc_doppl_1st = self.pilot_loc_doppl_1st;
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
            is_in = self.isInArea(1, pos_doppl, varargin{:});
        end
        function is_in = isInAreaCE(self, pos_doppl, varargin)
            is_in = self.isInArea(2, pos_doppl, varargin{:});
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
        insert pilots
        @pilots_pow: pilot power
        %}
        function insertP(self, pilots_pow)
            % input check
            self.pilots_len = length(self.pilots);
            % input check - pilot numbers
            if self.pilots_num_delay < 0 || self.pilots_num_doppl < 0
                error("Pilot number on each axis must be non-negative.");
            end
            % input check - pilots
            if ~isempty(self.pilots)
                if self.pilots_len ~= self.pilots_num_delay*self.pilots_num_doppl
                    error("The manual pilot input do not have the required number.");
                elseif self.pilots_len > self.nTimeslotNum*self.nSubcarNum
                    error("The manual pilot input overflows (over the OTFS frame size).");
                end
            end
            % input check - pilot power
            if self.pilots_len == 0 && isnan(pilots_pow)
                error("The pilots (linear) power is required while no manual pilot input.");
            end
            % initiate pilots if empty
            if isempty(self.pilots)
                self.pilots_len = self.pilots_num_delay*self.pilots_num_doppl;
                self.pilots = sqrt(pilots_pow/2)*(1+1j)*ones(self.pilots_len, 1);
            end
            % allocate pilots
            if self.pilots_len ~= 0
                % allocate the coordinate for 1st pilot
                switch self.pilot_loc_type
                    case self.PILOT_LOC_CENTER
                        self.pilot_loc_delay_1st = floor((self.nSubcarNum - self.pilots_num_delay)/2) + 1;
                        self.pilot_loc_doppl_1st = floor((self.nTimeslotNum - self.pilots_num_doppl)/2) + 1;
                    case self.PILOT_LOC_ZP 
                        self.pilot_loc_delay_1st = floor((self.nSubcarNum - self.pilots_num_delay)/2);
                        self.pilot_loc_doppl_1st = self.nTimeslotNum - self.zp_len + floor((self.zp_len - self.pilots_num_doppl)/2);
                end
                % allocate pilots
                self.content(self.pilot_loc_doppl_1st:self.pilot_loc_doppl_1st+self.pilots_num_doppl-1, self.pilot_loc_delay_1st:self.pilot_loc_delay_1st+self.pilots_num_delay-1) = transpose(reshape(self.pilots, self.pilots_num_delay, self.pilots_num_doppl));
            end
        end

        %{
        insert guards
        %}
        function insertG(self, guard_delay_full, guard_doppl_full, guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos)
            % input check - full guard
            if guard_delay_full
                % overwrite guard settings
                guard_delay_num_neg = floor((self.nSubcarNum - self.pilots_num_delay)/2);
                guard_delay_num_pos = self.nSubcarNum - self.pilots_num_delay - guard_delay_num_neg;
            else
                if guard_delay_num_neg < 0 || guard_delay_num_pos < 0
                    error("Guard number along the delay axis must be non-negative.");
                end
                if floor(guard_delay_num_neg) ~= guard_delay_num_neg || floor(guard_delay_num_pos) ~= guard_delay_num_pos
                    error("Guard number along the delay/Doppler axis must be integers.");
                end
            end
            if guard_doppl_full
                guard_doppl_num_neg = floor((self.nTimeslotNum - self.pilots_num_doppl)/2);
                guard_doppl_num_pos = self.nTimeslotNum - self.pilots_num_doppl - guard_doppl_num_neg;
            else
                if guard_doppl_num_neg < 0 || guard_doppl_num_pos < 0
                    error("Guard number along the Doppler axis must be non-negative.");
                end
                if floor(guard_doppl_num_neg) ~= guard_doppl_num_neg || floor(guard_doppl_num_pos) ~= guard_doppl_num_pos
                    error("Guard number along the Doppler axis must be integers.");
                end
            end
            % input check - overflow (guards)
            if self.pilot_loc_delay_1st - guard_delay_num_neg <= 0
                error("The guard (neg) on delay axis overflows.");
            end
            if (self.pilot_loc_delay_1st+self.pilots_num_delay-1) + guard_delay_num_pos > self.nSubcarNum
                error("The guard (pos) on delay axis overflows.");
            end
            if self.pilot_loc_doppl_1st - guard_doppl_num_neg <= 0
                error("The guard (neg) on Doppler axis overflows.");
            end
            if (self.pilot_loc_doppl_1st+self.pilots_num_doppl-1) + guard_doppl_num_pos > self.nTimeslotNum
                error("The guard (pos) on Doppler axis overflows.");
            end
            
            % calculate area
            if self.pilots_len > 0
                % calculate PG area
                self.pg_num = (self.pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(self.pilots_num_doppl+guard_doppl_num_neg+guard_doppl_num_pos);
                self.pg_delay_beg = self.pilot_loc_delay_1st - guard_delay_num_neg;
                self.pg_delay_end = self.pilot_loc_delay_1st + self.pilots_num_delay - 1 + guard_delay_num_pos;
                self.pg_doppl_beg = self.pilot_loc_doppl_1st - guard_doppl_num_neg;
                self.pg_doppl_end = self.pilot_loc_doppl_1st + self.pilots_num_doppl - 1 + guard_doppl_num_pos;
                % calulate channel estimate area
                if guard_delay_full
                    self.ce_delay_beg = 1;
                    self.ce_delay_end = self.nSubcarNum;
                else
                    self.ce_delay_beg = self.pg_delay_beg + guard_delay_num_neg;
                    self.ce_delay_end = self.pg_delay_end;
                end
                if guard_doppl_full
                    self.ce_doppl_beg = 1;
                    self.ce_doppl_end = self.nTimeslotNum;
                else
                    self.ce_doppl_beg = self.pg_doppl_beg + floor(guard_doppl_num_neg/2);
                    self.ce_doppl_end = self.pg_doppl_end - floor(guard_doppl_num_pos/2);
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
            data_num = self.nTimeslotNum*self.nSubcarNum - self.pg_num - self.zp_len*self.nTimeslotNum;
            % input check - symbols
            if ~isvector(symbols)
                error("The transmission symbol must be a vector");
            elseif length(symbols) ~= data_num
                error("The transmission symbol number must be %d", data_num);
            end
            
            % modulate
            % reshape(rowwise) to [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
            symbols = symbols(:);
            if data_num == self.nSubcarNum*self.nTimeslotNum
                self.content = transpose(reshape(symbols, self.nSubcarNum, self.nTimeslotNum));
            else
                symbols_id = 1;
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if ~self.isInAreaPG(doppl_id, delay_id)
                            self.content(doppl_id, delay_id) = symbols(symbols_id);
                            symbols_id = symbols_id + 1;
                        end
                    end
                end
                assert(symbols_id - 1 == data_num);
            end
        end

        %{
        decide whether the current location is in CE area
        @tag:           1->PG, 2->CE
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
            if tag == 1
                is_in = delay_pos >= self.pg_delay_beg && delay_pos <= self.pg_delay_end && pos_doppl >= self.pg_doppl_beg && pos_doppl <= self.pg_doppl_end;
            elseif tag == 2
                is_in = delay_pos >= self.ce_delay_beg && delay_pos <= self.ce_delay_end && pos_doppl >= self.ce_doppl_beg && pos_doppl <= self.ce_doppl_end;
            end
        end

        %{
        % estimated the channel using a sngle pilot
        % @threshold: the threshold to detect a path
        %}
        function [his, lis, kis] = estimateChannel4SingPilot(self, threshold)
            % input check
            if self.pilots_len ~= 1
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
                        li = delay_id - self.pilot_loc_delay_1st;
                        ki = doppl_id - self.pilot_loc_doppl_1st;
                        if self.pulse_type == self.PULSE_IDEAL
                            pss_beta = exp(-2j*pi*li*ki/self.nSubcarNum/self.nTimeslotNum);
                        elseif self.pulse_type == self.PULSE_RECTA
                            pss_beta = exp(2j*pi*(self.pilot_loc_delay_1st - 1)*ki/self.nSubcarNum/self.nTimeslotNum);
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