# header
BEGIN {
    Iteration=0
    resetCounters()
}

# reset counters used for variable postfix
function resetCounters() {
    clockTimeCnt=0
    executionTimeCnt=0
    phiCnt=0
    phiFinalResCnt=0
    phiItersCnt=0
    # Reset counters for general Solving for extraction
    for (varName in subIter)
    {
        subIter[varName]=0
    }
}

# Extract value after columnSel
function extract(inLine,columnSel,outVar,a,b)
{
    a=index(inLine, columnSel)
    b=length(columnSel)
    split(substr(inLine, a+b),outVar)
    gsub("[,:]","",outVar[1])
}

# Iteration separator (increments 'Iteration')
/^[ \t]*Time = / {
    Iteration++
    resetCounters()
}

# Time extraction (sets 'Time')
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    Time=val[1]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

    # Extraction of any solved for variable
    /Solving for/ {
        extract($0, "Solving for ", varNameVal)

        varName=varNameVal[1]
        file=varName "_" subIter[varName]++
        file="./logs/" file
        extract($0, "Initial residual = ", val)
        print Time "\t" val[1] > file

        varName=varNameVal[1] "FinalRes"
        file=varName "_" subIter[varName]++
        file="./logs/" file
        extract($0, "Final residual = ", val)
        print Time "\t" val[1] > file

        varName=varNameVal[1] "Iters"
        file=varName "_" subIter[varName]++
        file="./logs/" file
        extract($0, "No Iterations ", val)
        print Time "\t" val[1] > file
    }

# Extraction of clockTime
/ClockTime = / {
    extract($0, "ClockTime =", val)
    file="./logs/clockTime_" clockTimeCnt
    print Time "	" val[1] > file
    clockTimeCnt++
}

# Extraction of executionTime
/ExecutionTime = / {
    extract($0, "ExecutionTime = ", val)
    file="./logs/executionTime_" executionTimeCnt
    print Time "	" val[1] > file
    executionTimeCnt++
}

