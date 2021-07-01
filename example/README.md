# VT-stability

Generating stability diagram

```console
$ # from root dir of repo
$ cd example
$ julia --project=.. c1c5_vtstability_print.jl > /tmp/temp.tsv
$ gnuplot -p -e 'file="/tmp/temp.tsv"' -c plot/vtstability-conc-temp-range.gnuplot
```
