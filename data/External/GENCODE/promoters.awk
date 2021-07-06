BEGIN{
    FS=OFS="\t";
}
{
    if ($4 == "+") {
        print $1, max(0, $2 - 1500), $2 + 500, $4, $5, $6
    } else {
        print $1, max(0, $2 - 500), $2 + 1500, $4, $5, $6
    }
}
function max(a, b) {
    return a > b ? a: b
}