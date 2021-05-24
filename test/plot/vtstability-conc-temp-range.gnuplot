# usage
#   gnuplot -p -e 'file = "<path_to_vtstab_crash_test.tsv>"' -c <path_to_this_script>
# This scripts works only with binary mixtures, because it's hardcoded.

set datafile separator tab

# colums
temperature = 1
conc = 2
error = 5
isstable = 6

# functions
is_error(x) = x eq "true"
is_stable(x) = x eq "true"

set xlabel "overall concentration (mol m^-^3)"
set ylabel "temperature (К)"

set key outside box

plot file u 2:( ! is_error(stringcolumn(error)) && is_stable(stringcolumn(isstable)) ? column(temperature) : 1/0) ls 9 w p t "stable single phase state",\
     file u 2:( ! is_error(stringcolumn(error)) && ! is_stable(stringcolumn(isstable)) ? column(temperature) : 1/0) ls 10 w p t "unstable single phase state",\
     file u 2:( is_error(stringcolumn(error)) ? column(temperature) : 1/0) ls 7 w p t "solver's error"
