clear 
use Raw_Data_ECO418 
log using log_Part1_ECO418

/* ESTIMATION */
keep if cop == 0
keep if ap1 == "JFK" | ap1 == "EWR" | ap1 == "LGA"


sort ap2 time
by ap2 time: egen market_population = sum(pax)
gen market_share = pax/market_population


sort ap2 time ap1
by ap2 time ap1: egen airport_population = sum(pax)
by ap2 time ap1: egen number_airlines = count(cr1)
gen airport_share = pax/airport_population
gen log_airport_share = log(airport_share)


keep if cr1 == "DL" | cr1 == "UA" | cr1 == "AA"

sort ap2 time
by ap2 time: egen inside_market_share = sum(market_share)
gen outside_market_share = 1 - inside_market_share
gen log_s_minus_log_0 = log(market_share) - log(outside_market_share)

gen expense = expense_per_asm*nsdst
keep if expense <.

gen ap1_ap2_cr1 = ap1+"-"+ap2+"-"+cr1
tabulate ap1_ap2_cr1, generate(ap1_ap2_cr1_id)

ivreg log_s_minus_log_0 (avprc log_airport_share = expense number_airlines) ///
		ap1_ap2_cr1_id*, noc
predict predicted 
global alpha = _b[avprc]
global sigma = _b[log_airport_share]
gen fixed_effect = predicted - $alpha *avprc - $sigma *log_airport_share


/* INSTRUMENT RELEVANCE */		
reg avprc expense number_airlines ap1_ap2_cr1_id*, noc
reg log_airport_share expense number_airlines ap1_ap2_cr1_id*, noc
cor avprc expense
cor log_airport_share number_airlines

/* CLOSE */
export delimited intermediate_data_for_Part2.csv
log close