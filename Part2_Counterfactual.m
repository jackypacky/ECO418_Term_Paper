diary log_Part2_ECO418

global alpha sigma
alpha = -0.010995;
sigma = 0.3127083;
summary = [];
table = readtable("intermediate_data_for_Part2.csv");

table.time = strcat(string(table.yr), "-", string(table.qtr));
destinations = unique(table.ap2);
no_destinations = height(destinations);
for destination_index = 1:no_destinations
    destination = destinations(destination_index);
    destination_index + " of " + no_destinations
    carriers_route_table_times = table(ismember(table.ap2, destination), :);
    
    times = unique(carriers_route_table_times.time);
    for time_index = 1:height(times)
        time = times(time_index);
        carriers_route_table = carriers_route_table_times(ismember(carriers_route_table_times.time, time), :);
        counterfactual_route = [];
        
        carriers = unique(carriers_route_table.cr1);
        for carrier_index = 1:height(carriers)
            carrier_name = carriers(carrier_index);
            route_table = carriers_route_table(ismember(carriers_route_table.cr1, carrier_name), :);
            
            for i = 1:height(route_table)
                ap_i = string(route_table.ap1(i));
                share_derivatives = [];
            
                s_i = route_table(ismember(route_table.ap1, ap_i), :).market_share;
                for j = 1:height(route_table)
                    ap_j = string(route_table.ap1(j));
                    if ap_i == ap_j
                        ap_s_i = route_table(ismember(route_table.ap1, ap_i), :).airport_share;
            
                        share_derivatives(end+1) = alpha*(1/(1-sigma))*s_i*( ...
                            1-sigma*ap_s_i - (1-sigma)*s_i);
                    else
                        s_j = route_table(ismember(route_table.ap1, ap_j), :).market_share;
            
                        share_derivatives(end+1) = -alpha*s_i*s_j;
                    end
                end
                route_table.("share_derivative_" + ap_i) = share_derivatives.';
            end
            
            f = @(x)estimate_MCs(x, route_table);
            opts1 =  optimset('display','off');
            MCs = lsqnonlin(f, repmat(1e9, height(route_table), 1)', zeros(height(route_table), 1), [],[],[],[],[],[], opts1);
            
            counterfactual_route = [counterfactual_route; [route_table.cr1 ...
                route_table.ap1 array2table(MCs.') array2table(route_table.fixed_effect, 'VariableNames', {'fixed_effect'}) ...
                array2table(route_table.market_population, 'VariableNames', {'market_population'}) ...
                array2table(route_table.market_share, 'VariableNames', {'market_share'}) ...
                array2table(route_table.avprc, 'VariableNames', {'avprc'}) ...
                array2table(route_table.time, 'VariableNames', {'time'}) ...
                route_table.ap2
                ]];
        end
        counterfactual_route.Properties.VariableNames = ["cr1", "ap1", "MC", "fixed_effect", ...
            "market_population", "actual_market_share", "actual_price", "time", "ap2"];
        
        g = @(x)estimate_prices(x, counterfactual_route, alpha);
        opts2 =  optimset('display','off');
        number_carriers = height(unique(counterfactual_route.cr1));
        markups = lsqnonlin(g, ones(number_carriers, 1), [], [],[],[],[],[],[], opts2);
        counterfactual_route = calculate_market_shares(markups, counterfactual_route, alpha);
        summary = [summary; counterfactual_route];
    end
end

summary.actual_pax = summary.market_population.*summary.actual_market_share;
summary.counterfactual_pax = summary.market_population.*summary.counterfactual_market_share;
"Actual Price " + sum(summary.actual_pax.*summary.actual_price)/sum(summary.actual_pax)
"Counterfactual Price " + sum(summary.counterfactual_pax.*summary.counterfactual_price)/sum(summary.counterfactual_pax)
summary_MC = summary(summary.MC < 1e8, :);
"MC " + sum(summary_MC.actual_pax.*summary_MC.MC)/sum(summary_MC.actual_pax)
airports = ["JFK" "EWR" "LGA"];
for airport_index = 1:3
    airport = airports(airport_index)
    airport_summary = summary(ismember(summary.ap1, airport), :);
    sum(airport_summary.actual_pax.*airport_summary.actual_price)/sum(airport_summary.actual_pax)
    sum(airport_summary.counterfactual_pax.*airport_summary.counterfactual_price)/sum(airport_summary.counterfactual_pax)
    airport_summary_MC = airport_summary(airport_summary.MC < 1e8, :);
    sum(airport_summary_MC.actual_pax.*airport_summary_MC.MC)/sum(airport_summary_MC.actual_pax)
end
airlines = ["AA" "DL" "UA"];
for airline_index = 1:3
    airline = airlines(airline_index)
    airline_summary = summary(ismember(summary.cr1, airline), :);
    "Actual Price " + sum(airline_summary.actual_pax.*airline_summary.actual_price)/sum(airline_summary.actual_pax)
    "Counterfactual Price " + sum(airline_summary.counterfactual_pax.*airline_summary.counterfactual_price)/sum(airline_summary.counterfactual_pax)
    airline_summary_MC = airline_summary(airline_summary.MC < 1e8, :);
    "MC " + sum(airline_summary_MC.actual_pax.*airline_summary_MC.MC)/sum(airline_summary_MC.actual_pax)
end

function F = estimate_MCs(x, route_table)
    optimalities = [];
    for i = 1:height(route_table)
        ap_i = string(route_table.ap1(i));
    
        s_i = route_table(ismember(route_table.ap1, ap_i), :).market_share;
        optimality = s_i;
    
        for j = 1:height(route_table)
            ap_j = string(route_table.ap1(j));
            P = route_table(ismember(route_table.ap1, ap_j), :).avprc;
            c = x(j);
            share_derivative = route_table(ismember(route_table.ap1, ap_j), :).("share_derivative_" + ap_i);
            optimality = optimality + (P-c)*share_derivative;
        end
    
        optimalities(end+1) = optimality;
    end 
    F = optimalities;
end

function G = estimate_prices(x, counterfactual_route, alpha)
counterfactual_route = calculate_market_shares(x, counterfactual_route, alpha);
carriers = unique(counterfactual_route.cr1);
optimalities = [];
for carrier_index = 1:height(carriers)
    carrier_name = string(carriers(carrier_index));
    m = x(carrier_index);
    optimality = m + 1/(alpha*(1- ...
        sum(counterfactual_route(ismember(counterfactual_route.cr1, carrier_name), "counterfactual_market_share")).counterfactual_market_share));
    optimalities(end+1) = optimality;
end
G = optimalities;
end 

function calculate_shares = calculate_market_shares(x, counterfactual_route, alpha)
markups = [];
carriers = unique(counterfactual_route.cr1);
for carrier_index = 1:height(carriers)
    carrier_name = carriers(carrier_index);
    m = x(carrier_index);
    markups = [markups repmat(m, sum(ismember(counterfactual_route.cr1, carrier_name)), 1)'];
end
counterfactual_route.m = markups';
counterfactual_route.counterfactual_price = counterfactual_route.MC + ...
    counterfactual_route.m;
counterfactual_route.EV = exp(alpha*(counterfactual_route.counterfactual_price) ...
    + counterfactual_route.fixed_effect);
total_EV = 1 + sum(counterfactual_route.EV);
counterfactual_route.counterfactual_market_share = counterfactual_route.EV / total_EV;
calculate_shares = counterfactual_route;
end

diary off
type log_Part2_ECO418