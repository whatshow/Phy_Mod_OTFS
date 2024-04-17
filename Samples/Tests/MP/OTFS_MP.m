function x_est = OTFS_MP(N,M,p,his,lis,kis,sigma_2,y,constel)  
    M_mod = length(constel);
    n_ite = 200;
    delta_fra = 0.6;

    conv_rate_prev = -0.1;

    % init all parameters
    % mean & variance (d, c), d=M*(l-1)+k, c=M*(l-li-1)+(k-ki)
    mu_dc = zeros(N*M, p);
    sigma2_dc = zeros(N*M, p);
    % probability (d, c), d=M*(l-1)+k, c=M*(l-1+li)+(k+ki)
    p_dc = ones(N*M, p, M_mod)*(1/M_mod);
    for ite=1:n_ite
        % we update the mean and variance for each d in mu[d, c]
        for l=1:M
            for k=1:N
                d = N*(l-1)+k;
                %TODO: jump if y[d] is in CE area
                mu_d = zeros(p,1);                              % the sum of mu[d, c] for a given d (must be initialised as 0)
                sigma2_d = zeros(p,1);                          % the sum of sigma2[d, c] for a given d (must be initialised as 0)
                % we consider all x[e] -> y[d]
                for p_id=1:p
                    hi = his(p_id);
                    li = lis(p_id);
                    ki = kis(p_id);
                    % calculate coordinates for x[e]
                    e_l = l - 1 - li + 1;
                    e_k = mod(k - 1 - ki, N) + 1;
                    e_h = hi*exp(2j*(pi/M)*(e_l-1)*ki/N);
                    if l-1 < li 
                        e_l = e_l + M;
                        e_h = e_h*exp(-2j*pi*(e_k-1)/N);
                    end
                    %TODO: jump if x[e] is in PG area
                    % for the current x[c], we consider all constellation points
                    for i2=1:1:M_mod
                        mu_d(p_id) = mu_d(p_id) + p_dc(d,p_id,i2) * constel(i2);
                        sigma2_d(p_id) = sigma2_d(p_id) + p_dc(d,p_id,i2) * abs(constel(i2))^2;
                    end
                    mu_d(p_id) = mu_d(p_id) * e_h;
                    sigma2_d(p_id) = sigma2_d(p_id) * abs(e_h)^2  - abs(mu_d(p_id))^2;
                end
                mu_d_sum = sum(mu_d);
                sigma2_d_sum = sum(sigma2_d)+(sigma_2);
                % remove x[c] from the sum
                for p_id=1:p
                    mu_dc(d,p_id) = mu_d_sum - mu_d(p_id);
                    sigma2_dc(d,p_id) = sigma2_d_sum - sigma2_d(p_id);
                end
            end
        end
        % Update probabilities for each c in p[d, c]
        p_cj = zeros(N*M, M_mod);                               % Pc(aj): the probability that x[c] for each constellation point
        for l=1:1:M
            for k=1:1:N
                c = N*(l-1)+k;
                %TODO: jump if x[c] in PG area
                pr_ecj_ln = zeros(p, M_mod);                    % ln(Pr(y[e]|x[c]=aj, H))
                % we consider all y[e] from x[c]
                e_ls = zeros(p,1);                              % y[e] coordinate l
                e_ks = zeros(p,1);                              % y[e] coordinate k
                for p_id=1:p
                    hi = his(p_id);
                    li = lis(p_id);
                    ki = kis(p_id);
                    % y[e] - calculate coordinates 
                    e_l = l - 1 + li + 1;
                    e_k = mod(k - 1 + ki, N) + 1;
                    if l + li <= M
                        e_h = hi*exp(2j*(pi/M)*(l-1)*ki/N);
                    else
                        e_l = e_l - M;
                        e_h = hi*exp(2j*(pi/M)*(l-1-M)*ki/N)*exp(-2j*pi*(k-1)/N);
                    end
                    % y[e] - record coordinates
                    e_ls(p_id) = e_l;
                    e_ks(p_id) = e_k;
                    % calculate ln(eta(e,c,k)): the probability that x[c] takes a constellation point(k) based on y[e]
                    eta_ec_ln = zeros(M_mod,1);                 % ln(eta(e,c,k))
                    for i2=1:1:M_mod
                        eta_ec_ln(i2) = -(abs(y(e_k, e_l)- mu_dc(N*(e_l-1)+e_k,p_id) - e_h*constel(i2))^2)/sigma2_dc(N*(e_l-1)+e_k,p_id);
                    end
                    eta_ec_ln = eta_ec_ln - max(eta_ec_ln);     % subtract the maximal exponenet of ln(eta(e,c,k)) to keep mathematical stability
                    % calculate the probability that p[d, c] for the given d
                    pr_ecj_ln(p_id, :) = eta_ec_ln - log(sum(exp(eta_ec_ln)));
                end
                % calculate sum(ln(Pr(y[e]|x[c]=aj, H))) for e in J(c)
                pr_cj_ln = zeros(1, M_mod);                     % sum(ln(Pr(y[e]|x[c]=aj, H)))
                for i2=1:1:M_mod
                    pr_cj_ln(i2) = sum(pr_ecj_ln(:,i2));
                end
                % calculate p_cj = exp(pr_cj_ln)
                p_cj_all = exp(pr_cj_ln - max(pr_cj_ln));       % stay mathematical stability
                p_cj(c, :) = p_cj_all/sum(p_cj_all);            % normalise the sum of probability
                % calculate p_dc
                for p_id=1:1:p
                    e_l = e_ls(p_id);
                    e_k = e_ks(p_id);
                    % remove ln(Pr(y[d]|x[c]=aj, H)) from sum(ln(Pr(y[e]|x[c]=aj, H)))
                    p_dc_val = pr_cj_ln - pr_ecj_ln(p_id,:);
                    p_dc_val = exp(p_dc_val - max(p_dc_val));
                    p_dc_val = p_dc_val/sum(p_dc_val);
                    p_dc(N*(e_l-1)+e_k,p_id,:) = p_dc_val*delta_fra + (1-delta_fra)*reshape(p_dc(N*(e_l-1)+e_k,p_id,:),1,M_mod);
                end

            end
        end
        % early stop
        conv_rate =  sum(max(p_cj,[],2)>0.99)/(N*M);
        if conv_rate==1
            sum_prob_fin = p_cj;
            break;
        elseif conv_rate > conv_rate_prev
            conv_rate_prev = conv_rate;
            sum_prob_fin = p_cj;
        elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
            break;
        end
    end
    x_est = zeros(N,M);
    for l=1:1:M
        for k=1:1:N
            [~,pos] = max(sum_prob_fin(N*(l-1)+k,:));
            x_est(k,l) = constel(pos);
        end
    end
end