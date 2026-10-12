import math

def checked_counts(counts):
    counts=tuple(counts)
    k=len(counts).bit_length()-1
    if k not in (3,4) or len(counts)!=2**k:
        raise ValueError('Expected 8 counts for 3 CpGs or 16 counts for 4 CpGs.')
    if any(not isinstance(n,(int,float)) or not math.isfinite(n) or n<0 or n!=int(n) for n in counts):
        raise ValueError('Counts must be finite nonnegative integers.')
    return counts,k

def metrics(counts):
    counts,k=checked_counts(counts);n=sum(counts)
    if n==0:return {'n_complete_reads':0,'beta':None,'entropy':None}
    beta=sum(i.bit_count()*v for i,v in enumerate(counts))/(k*n)
    h=-sum((v/n)*math.log2(v/n) for v in counts if v)/k
    return {'n_complete_reads':int(n),'beta':beta,'entropy':h}

def compare_counts(control,test):
    control,k=checked_counts(control);test,j=checked_counts(test)
    if k!=j:raise ValueError('Comparison requires identical CpG-unit size.')
    a=metrics(control);b=metrics(test)
    if not sum(control) or not sum(test):
        return {'delta_beta':None,'delta_entropy':None,'D_epi':None}
    p=[v/sum(control) for v in control];q=[v/sum(test) for v in test]
    js=0.
    for x,y in zip(p,q):
        m=(x+y)/2
        if x:js+=.5*x*math.log2(x/m)
        if y:js+=.5*y*math.log2(y/m)
    return {'delta_beta':b['beta']-a['beta'],'delta_entropy':b['entropy']-a['entropy'],'D_epi':math.sqrt(max(0.,min(1.,js)))}
