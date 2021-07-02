# usage
#   gnuplot -p -e 'file = "<path_to_vtstab_crash_test.tsv>"' -c <path_to_this_script>

set datafile separator tab

# columns
temperature = 1
conc = 2
converged = 3
singlephase = 4

is_converged(x) = x eq "true"
is_singlephase(x) = x eq "true"

plot file u conc:( is_converged(stringcolumn(converged)) && is_singlephase(stringcolumn(singlephase)) ? column(temperature) : 1/0 ) ls 9 w p t "1-phase",\
    file u conc:( is_converged(stringcolumn(converged)) && ! is_singlephase(stringcolumn(singlephase)) ? column(temperature) : 1/0 ) ls 10 w p t "2-phase",\
     file u conc:( ! is_converged(stringcolumn(converged)) ? column(temperature) : 1/0 ) ls 7 w p t "not converged"
