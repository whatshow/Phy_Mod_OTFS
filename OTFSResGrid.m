% resource grid
classdef OTFSResGrid < handle
    properties(Constant)
        % Pilot locations
        PILOT_LOC_CENTER = 10;                      % the pilot is put at the center of frame
        PILOT_LOC_DELAY_MOST_CENTER = 20;           % the pilot is put at the most delay area center
        PILOT_LOC_TYPES = [OTFS.PILOT_LOC_CENTER, OTFS.PILOT_LOC_DELAY_MOST_CENTER];
    end
    properties
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        X_DD
        % zero padding
        zp_len = 0;
        % pilot
        pilots = [];                                % pilots value
        pilots_num_delay = 0;                       % pilots number along the delay(Doppler) axis
        pilots_num_doppl = 0;                       % pilots number along the Doppler(time) axis
        pilot_loc_delay_1st = 0;                    % 1st (lowest) pilot location in delay axis
        pilot_loc_doppl_1st = 0;                    % 1st (lowest) pilot location in Doppler axis
        % pilot location
        pilot_loc_type = OTFS.PILOT_LOC_CENTER;
        % channel estimation area in X_DD
        ce_xdd_num = 0;
        ce_xdd_delay_beg = NaN;
        ce_xdd_delay_end = NaN;
        ce_xdd_doppl_beg = NaN;
        ce_xdd_doppl_end = NaN;
        % channel estimation area in Y_DD
        ce_ydd_delay_beg = NaN;
        ce_ydd_delay_end = NaN;
        ce_ydd_doppl_beg = NaN;
        ce_ydd_doppl_end = NaN;
    end
    methods
        %{
        init the resource grid
        %}
        function self = OTFSResGrid(nSubcarNum, nTimeslotNum, varargin)
            % optional inputs
            ip = inputParser;
            addParameter(ip, "pilot_loc_type", self.pilot_loc_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, self.PILOT_LOC_TYPES));
            addParameter(ip, "zp_len", self.zp_len, @(x) isscalar(x)&&isnumeric(x));
            ip.KeepUnmatched = true;     % Allow unmatched cases
            ip.CaseSensitive = false;    % Allow capital or small characters
            parse(ip, varargin{:});
            % take inputs
            % nSubcarNum  
            if ~isscalar(nSubcarNum) || nSubcarNum ~= floor(nSubcarNum)
                error("The number of subcarriers must be an integer scalar.");
            else
                self.nSubcarNum = nSubcarNum;
            end
            % nSubcarNum
            if ~isscalar(nTimeslotNum) || nTimeslotNum ~= floor(nTimeslotNum)
                error("The number of timeslots must be an integer scalar.");
            else
                self.nTimeslotNum = nTimeslotNum;
            end
            % pilot location
            self.pilot_loc_type = ip.Results.pilot_loc_type;
            % zero padding
            self.zp_len = ip.Results.zp_len;
            if self.zp_len < 0 || self.zp_len > self.nSubcarNum
                error("Zero Padding length cannot be negative or over subcarrier number.");
            end
            % init X_DD
            self.X_DD = zeros(self.nTimeslotNum, self.nSubcarNum);
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
        @force:                     whether we force to insert pilots when there is no enough space
        %}
        function map(self, symbols, varargin)
            % optional inputs - register
            ip = inputParser;
            addParameter(ip, 'pilots', self.pilots, @(x) isvector(x)&&isnumeric(x));
            addParameter(ip, 'pilots_pow', NaN, @(x) isscalar(x)&&isnumeric(x));
            addParameter(ip, 'pilots_num_delay', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(ip, 'pilots_num_doppl', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(ip, 'guard_delay_full', false, @(x) isscalar(x)&&islogical(x));
            addParameter(ip, 'guard_doppl_full', false, @(x) isscalar(x)&&islogical(x));
            addParameter(ip, 'guard_delay_num_neg', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(ip, 'guard_delay_num_pos', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(ip, 'guard_doppl_num_neg', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(ip, 'guard_doppl_num_pos', 0, @(x) isscalar(x)&&isnumeric(x));
            ip.KeepUnmatched = true;     % Allow unmatched cases
            ip.CaseSensitive = false;    % Allow capital or small characters
            parse(ip, varargin{:});
            % optional inputs - assign
            self.pilots              = ip.Results.pilots;
            self.pilots_num_delay    = ip.Results.pilots_num_delay;
            self.pilots_num_doppl    = ip.Results.pilots_num_doppl;
            pilots_pow          = ip.Results.pilots_pow;
            guard_delay_full    = ip.Results.guard_delay_full;
            guard_doppl_full    = ip.Results.guard_doppl_full;
            guard_delay_num_neg = ip.Results.guard_delay_num_neg;
            guard_delay_num_pos = ip.Results.guard_delay_num_pos;
            guard_doppl_num_neg = ip.Results.guard_doppl_num_neg;
            guard_doppl_num_pos = ip.Results.guard_doppl_num_pos;

            % insert pilots and guards
            self.insertPG(pilots_pow, guard_delay_full, guard_doppl_full, guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos);
            % insert data
            self.insertDA(symbols);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %{
        get X_DD
        %}
        function X_DD = getXDD(self)
            X_DD = self.X_DD;
        end

        %{
        check whether the pilots & guards are inserted or not
        %}
        function have = havePG(self)
            % check whether pilots is assigned or not
            have = ~isempty(self.pilots);
            % check whether X_DD_invalid area is calculated
            have = have && ~isnan(self.ce_xdd_delay_beg) && ~isnan(self.ce_xdd_delay_end) && ~isnan(self.ce_xdd_doppl_beg) && ~isnan(self.ce_xdd_doppl_end);
            % check whether CE area is calulated
            have = have && ~isnan(self.ce_ydd_delay_beg) && ~isnan(self.ce_ydd_delay_end) && ~isnan(self.ce_ydd_doppl_beg) && ~isnan(self.ce_ydd_doppl_end);
        end

        %{
        check whether the current position is in channel estimation area
        @ce_area_tag:   notify whether CE area to check
        @pos_delay:     the position on the delay axis for matrix or th position on the Doppler-delay axis for the vector
        @pos_doppl:     the position on the Doppler axis. Not given means the position is for a Doppler-delay vector
        %}
        function is_in = isInCEAreaXDD(self, delay_pos, varargin)
            is_in = self.isInCEArea(1, delay_pos, varargin{:});
        end
        function is_in = isInCEAreaYDD(self, delay_pos, varargin)
            is_in = self.isInCEArea(2, delay_pos, varargin{:});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access=private)
        %{
        insert pilots and guards
        %}
        function insertPG(self, pilots_pow, guard_delay_full, guard_doppl_full, guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos)
            % input check
            pilots_len = length(self.pilots);
            % input check - pilot numbers
            if self.pilots_num_delay < 0 || self.pilots_num_doppl < 0
                error("Pilot number on each axis must be non-negative.");
            end
            % input check - pilots
            if ~isempty(self.pilots)
                if pilots_len ~= self.pilots_num_delay*self.pilots_num_doppl
                    error("The manual pilot input do not have the required number.");
                elseif pilots_len > self.nTimeslotNum*self.nSubcarNum
                    error("The manual pilot input overflows (over the OTFS frame size).");
                end
            end
            % input check - pilot power
            if pilots_len == 0 && isnan(pilots_pow)
                error("The pilots (linear) power is required while no manual pilot input.");
            end
            % input check - guard
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
            % input check - overflow (pilots + guards)
            if self.pilots_num_delay + guard_delay_num_neg + guard_doppl_num_pos > self.nSubcarNum
                error("Overflow on the delay axis (pilots + guards is over subcarrier number).");
            end
            if self.pilots_num_doppl + guard_doppl_num_neg + guard_doppl_num_pos > self.nTimeslotNum
                error("Overflow on the Doppler axis (pilots + guards is over timeslot number).");
            end

            % initiate pilots if empty
            if isempty(self.pilots)
                pilots_len = self.pilots_num_delay*self.pilots_num_doppl;
                self.pilots = sqrt(pilots_pow/2)*(1+1j)*ones(pilots_len, 1);
            end
            % allocate pilots
            if pilots_len == 0
                % no pilots no operation
                self.ce_xdd_num = 0;
            else
                % some pilots
                % calulate the data number for two axises
                data_delay_num = (self.nSubcarNum - self.pilots_num_delay - guard_delay_num_neg - guard_delay_num_pos);
                data_doppler_num = self.nTimeslotNum - self.pilots_num_doppl - guard_doppl_num_neg - guard_doppl_num_pos;
                % calculate the pilot start point shift due to asymmetric guards
                pilot_shift_delay_pos = guard_delay_num_pos - guard_delay_num_neg; 
                pilot_shift_doppler_pos = guard_doppl_num_pos - guard_doppl_num_neg;
                % locate the 1st coordinates of the pilots (following the positive direction of axises)
                self.pilot_loc_delay_1st = 0;
                self.pilot_loc_doppl_1st = 0;
                switch self.pilot_loc_type
                    case self.PILOT_LOC_CENTER
                        self.pilot_loc_delay_1st = floor(data_delay_num/2) + guard_delay_num_neg + pilot_shift_delay_pos + 1;
                        self.pilot_loc_doppl_1st = floor(data_doppler_num/2) + guard_doppl_num_neg + pilot_shift_doppler_pos + 1;
                    case self.PILOT_LOC_DELAY_MOST_CENTER 
                        self.pilot_loc_delay_1st = floor(data_delay_num/2) + guard_delay_num_neg + pilot_shift_delay_pos + 1;
                        self.pilot_loc_doppl_1st = data_doppler_num + guard_doppl_num_neg + 1;
                end
                % allocate pilots
                self.X_DD(self.pilot_loc_doppl_1st:self.pilot_loc_doppl_1st+self.pilots_num_doppl-1, self.pilot_loc_delay_1st:self.pilot_loc_delay_1st+self.pilots_num_delay-1) = transpose(reshape(self.pilots, self.pilots_num_delay, self.pilots_num_doppl));
                % calculate the invalid area in X_DD
                self.ce_xdd_num = (self.pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(self.pilots_num_doppl+guard_doppl_num_neg+guard_doppl_num_pos);
                self.ce_xdd_delay_beg = self.pilot_loc_delay_1st - guard_delay_num_neg;
                self.ce_xdd_delay_end = self.pilot_loc_delay_1st + self.pilots_num_delay - 1 + guard_delay_num_pos;
                self.ce_xdd_doppl_beg = self.pilot_loc_doppl_1st - guard_doppl_num_neg;
                self.ce_xdd_doppl_end = self.pilot_loc_doppl_1st + self.pilots_num_doppl - 1 + guard_doppl_num_pos;
            end
            % calulate channel estimate area
            if pilots_len > 0
                if guard_delay_full
                    self.ce_ydd_delay_beg = 1 + guard_doppl_num_neg;
                    self.ce_ydd_delay_end = self.nSubcarNum;
                else
                    self.ce_ydd_delay_beg = self.ce_xdd_delay_beg + guard_delay_num_neg;
                    self.ce_ydd_delay_end = self.ce_xdd_delay_end;
                end
                if guard_doppl_full
                    self.ce_ydd_doppl_beg = 1 + floor(guard_doppl_num_neg/2);
                    self.ce_ydd_doppl_end = self.nTimeslotNum - floor(guard_doppl_num_pos/2);
                else
                    self.ce_ydd_doppl_beg = self.ce_xdd_doppl_beg + floor(guard_doppl_num_neg/2);
                    self.ce_ydd_doppl_end = self.ce_xdd_doppl_end - floor(guard_doppl_num_pos/2);
                end
            end
        end
        
        %{
        insert data
        @symbols: symbols to map a vector
        %}
        function insertDA(self, symbols)
            % input check
            data_num = self.nTimeslotNum*self.nSubcarNum - self.ce_xdd_num - self.zp_len*self.nTimeslotNum;
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
                self.X_DD = transpose(reshape(symbols, self.nSubcarNum, self.nTimeslotNum));
            else
                symbols_id = 1;
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if ~self.isInCEAreaXDD(delay_id, doppl_id)
                            self.X_DD(doppl_id, delay_id) = symbols(symbols_id);
                            symbols_id = symbols_id + 1;
                        end
                    end
                end
                assert(symbols_id - 1 == data_num);
            end
        end

        %{
        decide whether the current location is in CE area
        @ce_area_tag: 1->X_DD, 2->Y_DD
        @delay_pos:     the delay position
        @pos_doppl:     the position on the Doppler axis. Not given means the position is for a Doppler-delay vector
        %}
        function is_in = isInCEArea(self, ce_area_tag, delay_pos, varargin)
            % input check
             if isempty(varargin)
                if delay_pos <= 0 || delay_pos > self.nSubcarNum*self.nTimeslotNum
                    error("The vector position is out of the OTFS size.");
                end
            else
                if delay_pos <= 0 || delay_pos > self.nSubcarNum
                    error("The delay position is out of the subcarrier number.");
                end
                if varargin{1} <= 0 || varargin{1} > self.nTimeslotNum
                    error("The Doppler position is out of the timeslot number.");
                end
            end
            
            % recalculate the position
            if isempty(varargin)
                pos_doppl = delay_pos/self.nSubcarNum;
                if pos_doppl == floor(pos_doppl)
                    pos_doppl = floor(pos_doppl);
                else
                    pos_doppl = floor(pos_doppl) + 1;
                end
                delay_pos = delay_pos - (pos_doppl-1)*self.nSubcarNum;
            else
                pos_doppl = varargin{1};
            end
            % decide
            switch ce_area_tag
                case 1
                    is_in = delay_pos >= self.ce_xdd_delay_beg && delay_pos <= self.ce_xdd_delay_end && pos_doppl >= self.ce_xdd_doppl_beg && pos_doppl <= self.ce_xdd_doppl_end;
                case 2
                    is_in = delay_pos >= self.ce_ydd_delay_beg && delay_pos <= self.ce_ydd_delay_end && pos_doppl >= self.ce_ydd_doppl_beg && pos_doppl <= self.ce_ydd_doppl_end;
            end
        end
    end
end