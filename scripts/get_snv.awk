#!/bin/awk -f

{
    if ($0~"#" || length($4)==1 && length($5)==1) {
        print
    }
}