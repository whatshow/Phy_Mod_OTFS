function x_est = OTFS_MP_Embed(N,M,p,his,lis,kis,sigma_2,y,constel)  
    M_mod = length(constel);
    yv = y(:);
    n_ite = 200;
    delta_fra = 0.6;

    mean_int = zeros(N*M,p);
    var_int = zeros(N*M,p);
    p_map = ones(N*M,p,M_mod)*(1/M_mod);

    conv_rate_prev = -0.1;

    % init all parameters
    % mean & variance (d, c), d=M*(l-1)+k, c=M*(l-li-1)+(k-ki)
    mu_dc = zeros(N*M, p);
    sigma2_dc = zeros(N*M, p);
    % probability (d, c), d=M*(l-1)+k, c=M*(l-li-1)+(k-ki)
    p_dc = ones(N*M,p,M_mod)*(1/M_mod);

    for ite=1:n_ite
        % we update the mean and variance for each d in μ[d, c]
        %% Update mean and var
        for l=1:M
            for k=1:N
                mean_int_hat = zeros(p,1);
                var_int_hat = zeros(p,1);
                for p_id=1:p
                    m = l-1-lis(p_id)+1;
                    add_term = exp(1i*2*(pi/M)*(m-1)*(kis(p_id)/N));
                    add_term1 = 1;
                    if l-1<lis(p_id)
                        n = mod(k-1-kis(p_id),N) + 1;
                        add_term1 = exp(-1i*2*pi*((n-1)/N));
                    end
                    new_chan = add_term * (add_term1) * his(p_id);

                    for i2=1:1:M_mod
                        mean_int_hat(p_id) = mean_int_hat(p_id) + p_map(N*(l-1)+k,p_id,i2) * constel(i2);
                        var_int_hat(p_id) = var_int_hat(p_id) + p_map(N*(l-1)+k,p_id,i2) * abs(constel(i2))^2;
                    end
                    mean_int_hat(p_id) = mean_int_hat(p_id) * new_chan;
                    var_int_hat(p_id) = var_int_hat(p_id) * abs(new_chan)^2  - abs(mean_int_hat(p_id))^2;
                end

                mean_int_sum = sum(mean_int_hat);
                var_int_sum = sum(var_int_hat)+(sigma_2);

                for p_id=1:p
                    mean_int(N*(l-1)+k,p_id) = mean_int_sum - mean_int_hat(p_id);
                    var_int(N*(l-1)+k,p_id) = var_int_sum - var_int_hat(p_id);
                end

            end
        end
        %% Update probabilities
        sum_prob_comp = zeros(N*M,M_mod);
        dum_eff_ele1 = zeros(p,1);
        dum_eff_ele2 = zeros(p,1);
        for ele1=1:1:M
            for ele2=1:1:N
                dum_sum_prob = zeros(M_mod,1);
                log_te_var = zeros(p,M_mod);
                for p_id=1:p

                    if ele1+lis(p_id)<=M
                        eff_ele1 = ele1 + lis(p_id);
                        add_term = exp(1i*2*(pi/M)*(ele1-1)*(kis(p_id)/N));
                        int_flag = 0;
                    else
                        eff_ele1 = ele1 + lis(p_id)- M;
                        add_term = exp(1i*2*(pi/M)*(ele1-1-M)*(kis(p_id)/N));
                        int_flag = 1;
                    end
                    add_term1 = 1;
                    if int_flag==1
                        add_term1 = exp(-1i*2*pi*((ele2-1)/N));
                    end
                    eff_ele2 = mod(ele2-1+kis(p_id),N) + 1;
                    new_chan = add_term * add_term1 * his(p_id);

                    dum_eff_ele1(p_id) = eff_ele1;
                    dum_eff_ele2(p_id) = eff_ele2;
                    for i2=1:1:M_mod
                        dum_sum_prob(i2) = abs(yv(N*(eff_ele1-1)+eff_ele2)- mean_int(N*(eff_ele1-1)+eff_ele2,p_id) - new_chan * constel(i2))^2;
                        dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(N*(eff_ele1-1)+eff_ele2,p_id));
                    end
                    dum_sum = dum_sum_prob - max(dum_sum_prob);
                    dum1 = sum(exp(dum_sum));
                    log_te_var(p_id,:) = dum_sum - log(dum1);
                end
                for i2=1:1:M_mod
                    ln_qi(i2) = sum(log_te_var(:,i2));
                end
                dum_sum = exp(ln_qi - max(ln_qi));
                dum1 = sum(dum_sum);
                sum_prob_comp(N*(ele1-1)+ele2,:) = dum_sum/dum1;
                for p_id=1:1:p
                    eff_ele1 = dum_eff_ele1(p_id);
                    eff_ele2 = dum_eff_ele2(p_id);

                    dum_sum = log_te_var(p_id,:);
                    ln_qi_loc = ln_qi - dum_sum;
                    dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                    dum1 = sum(dum_sum);
                    p_map(N*(eff_ele1-1)+eff_ele2,p_id,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(N*(eff_ele1-1)+eff_ele2,p_id,:),1,M_mod);
                end

            end
        end
        conv_rate =  sum(max(sum_prob_comp,[],2)>0.99)/(N*M);
        if conv_rate==1
            sum_prob_fin = sum_prob_comp;
            break;
        elseif conv_rate > conv_rate_prev
            conv_rate_prev = conv_rate;
            sum_prob_fin = sum_prob_comp;
        elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
            break;
        end
    end
    x_est = zeros(N,M);
    for ele1=1:1:M
        for ele2=1:1:N
            [~,pos] = max(sum_prob_fin(N*(ele1-1)+ele2,:));
            x_est(ele2,ele1) = constel(pos);
        end
    end
end