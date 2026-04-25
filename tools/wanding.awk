function rev(string, x,i) {
		for(i=length;i!=0;i--)
				x=(x substr(string,i,1))
		return x
}

function count(string, c,    i, a, cnt) {
    cnt = 0
    for (i=0;i<n=split(string,a,"");i++)
        if (a[i]==c)
            cnt++
    return cnt;
}


function round(x,   ival, aval, fraction) {
   ival = int(x)    # integer part, int() truncates

   # see if fractional part
   if (ival == x)   # no fraction
      return ival   # ensure no decimals

   if (x < 0) {
      aval = -x     # absolute value
      ival = int(aval)
      fraction = aval - ival
      if (fraction >= .5)
         return int(x) - 1   # -2.5 --> -3
      else
         return int(x)       # -2.3 --> -2
   } else {
      fraction = x - ival
      if (fraction >= .5)
         return ival + 1
      else
         return ival
   }
}

function max(a,b) {
		if (a>b) 
				return a;
		else
				return b;
}

function min(a,b) {
		if (a<b) 
				return a;
		else
				return b;
}


function abs(value) {
    return (value<0?-value:value);
}

function minabs(a,b,  aa, ab) {
    aa = abs(a);
    ab = abs(b);
    if (aa < ab)
        return a;
    else
        return b;
}

function samIsRead1(flag) {
    return and(flag, 0x40);
}

function samIsRead2(flag) {
    return and(flag, 0x80);
}

function div(denom, divd) {
    return (divd == 0? "NA" : denom/divd);
}

function dnarev(string,        i, chars, dic, new_string) {
    dic["A"] = "T";
    dic["T"] = "A";
    dic["G"] = "C";
    dic["C"] = "G";
    dic["Y"] = "R";
    dic["R"] = "Y";
    dic["-"] = "-";
    dic["."] = ".";

    new_string = "";
    split(string, chars, "");
    for (i=length(chars);i>=0;i--) {
        if (chars[i] in dic)
            new_string = new_string dic[chars[i]];
        else
            new_string = new_string chars[i];
    }

    return new_string;
}

function aa3to1(string,    dic) {
    dic["Ala"]="A";
    dic["Arg"]="R";
    dic["Asn"]="N";
    dic["Asp"]="D";
    dic["Cys"]="C";
    dic["Glu"]="E";
    dic["Gln"]="Q";
    dic["Gly"]="G";
    dic["His"]="H";
    dic["Ile"]="I";
    dic["Leu"]="L";
    dic["Lys"]="K";
    dic["Met"]="M";
    dic["Phe"]="F";
    dic["Pro"]="P";
    dic["Ser"]="S";
    dic["Thr"]="T";
    dic["Trp"]="W";
    dic["Tyr"]="Y";
    dic["Val"]="V";
    dic["Sec"]="U";
    dic["Pyl"]="O";
    dic["X"]="*";
    return dic[string];
}

# join fields by specifying a range
# for example, from $3 to $10, one can use joinr(3,10)
function joinr(start, end,    sep, result, i) {

    if (sep == "")
        sep = "\t"
    result = ""
    for (i=start; i<=end; ++i) {
        if (i != start) {
            result = result sep;
        }
        result = result $i;
    }
    return result
}

function rep(s, n, sep,   i) {
    if (sep == "")
        sep = "\t"
    result = ""
    for (i=1; i<=n; ++i) {
        if (i != start) {
            result = result sep;
        }
        result = result s;
    }
    return result
}

function sumr(start, end,    sep, result, i) {

    if (sep == "")
        sep = "\t"
    result = 0
    for (i=start; i<=end; ++i) {
        if (i != start) {
            result = result sep;
        }
        result = result+$i;
    }
    return result
}

function minr(start,end,   first,result,i) {

    first = 1
    for (i=start; i<=end; ++i) {
        if ($i == ".") continue
        if (first || $i < result) {
            result = $i
            first = 0
        }
    }
    return result;
}

function maxr(start,end,   first,result,i) {

    first = 1
    for (i=start; i<=end; ++i) {
        if ($i == ".") continue
        if (first || $i > result) {
            result = $i
            first = 0
        }
    }
    return result;
}

function join(array, start, end, sep,  result, i)
{
    if (sep == "")
       sep = " "
    else if (sep == SUBSEP) # magic value
       sep = ""

    if (start == "")
        start = 1
    if (end == "")
        end = length(array)
    result = array[start]
    for (i = start + 1; i <= end; i++)
        result = result sep array[i]
    return result
}

function uniqarray(s, sep,       resultarr, result, i)
{
    split(s, sa, sep);
    split("", resultarr);
    for (i=1; i<=length(sa); i++)
        resultarr[sa[i]]=1;
    result=""
    for (i in resultarr) {
        if (result!="") result=result""sep
        result=result""i
    }
    return result
}


function joina(array, sep, start, end,    result, i)
{
    if (sep == "")
       sep = " "
    else if (sep == SUBSEP) # magic value
       sep = ""

    if (start == "")
        start = 1
    if (end == "")
        end = length(array)
    result = array[start]
    for (i = start + 1; i <= end; i++)
        result = result sep array[i]
    return result
}

# escape the string when the string is
# used as a pattern in regular expression
# i.e., \ needs to be replaced by \\
function escape(pat,   safe) {
  safe = pat
  gsub(/[][^$.*?+{}\\()|]/, "\\\\&", safe)
  return safe
}



function resplit(       a, l, i, j, b, k, BNF) # all are local variables
{
    ## l=split($0, a, /[\"\']/)
    l=split($0, a, /[\\"\\']/)
    BNF=0
    delete B
    for (i=1;i<=l;++i)
    {
        if (i % 2)
        {
            # print i,"1", a[i]
            k=split(a[i], b)
            for (j=1;j<=k;++j) {
                B[++BNF] = b[j]
                # print BNF, B[BNF];
            }
        }
        else
        {
            # print i,"0", a[i]
            B[++BNF] = a[i]
            # print BNF, B[BNF]
        }
    }
    # print "length", length(B)
    # print "\n"
}

function repeat(str, n, rep,    i )
{
    for(i=1; i<=n; i++)
        rep[i] = str;
    return
}

function repeats( str, n,    rep, i )
{
    for( ; i<n; i++ )
        rep = rep str   
    return rep
}

function labcol() {
  for(i=1;i<=NF;++i){$i=i"_"$i;}print $0;
}

function basename(file) {
    sub(".*/", "", file)
    return file
}

function basename2(file, a, n) {
  n = split(file, a, "/")
    return a[n]
}

function capFirst(s) {
	return toupper(substr(s,1,1))tolower(substr(s,2))
}


function FindAllMatches(str, regex, match_arr) {

    ftotal = 0;
    ini = RSTART;
    leng = RLENGTH;

    delete match_arr;

    while (match(str, regex) > 0) {
        match_arr[++ftotal] = substr(str, RSTART, RLENGTH)
        if (str == "") break
        str = substr(str, RSTART + (RLENGTH ? RLENGTH : 1))
    }

    RSTART = ini;
    RLENGTH = leng;
}
