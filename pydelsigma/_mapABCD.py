# pydelsigma/_mapABCD.py
from __future__ import division
import numpy as np

from ._utils import diagonal_indices

def mapABCD(ABCD, form='CRFB'):
    order = ABCD.shape[0] - 1
    odd = order % 2
    even = 1 - odd
    diagonal = diagonal_indices(ABCD[:order, :order])
    subdiag = diagonal_indices(ABCD[:order, :order], -1)
    supdiag = diagonal_indices([ABCD[:order, :order+1])(2 + odd-1):2:order] - 1
    if form in 'CRFB','CIFB','CRFBD':
        c = ABCD[subdiag]
        g = -ABCD[supdiag]
        if form == 'CRFB':
            dly = np.arange(1 + odd, order, 2)
            ABCD[dly, :] = ABCD[dly, :] - diag(c[dly - 1]) * ABCD[dly - 1,:]
        elif form == 'CRFBD':
            dly = np.arange(odd, order, 2)
            ABCD[dly, :] = ABCD[dly, :] + diag(g)*ABCD[dly + 1, :]
            if order > 2:
                coupl = np.arange(1 + even, order, 2)
                ABCD[coupl, :] = ABCD[coupl, :] - diag(c[coupl - 1]) * ABCD[coupl - 1, :]
        a = -ABCD[:order, order + 1].T
        b = ABCD[:, order].T
    elif 'CRFF' == form:
        c=np.array([- ABCD[0,(order + 2-1)],ABCD[(subdiag[(1-1):end - 1]-1)]]).reshape(1,-1)
        g=- ABCD[(supdiag-1)]
        if even:
            multg=range(1,(order+1),2)
            ABCD[(multg-1),:]=ABCD[(multg-1),:] + diag(g) * ABCD[(multg + 1-1),:]
        multc=range(3,(order+1),2)
        ABCD[(multc-1),:]=ABCD[(multc-1),:] - diag(c[(multc-1)]) * ABCD[(multc - 1-1),:]
        a(range(2,(order+1),2))=ABCD[(order + 1-1),(2-1):2:order]
        for i in range(2,(order+1),2):
            ABCD[(order + 1-1),:]=ABCD[(order + 1-1),:] - a[(i-1)] * ABCD[(i-1),:]
        a[(1-1):2:order]=ABCD[(order + 1-1),(1-1):2:order]
        b=ABCD[:,(order + 1-1)].T
    elif 'CRFFD' == form:
        #order=order - 1
        #odd=rem(order,2)
        #even=not  odd
        diagonal=diagonal[(1-1):order]
        subdiag=diagonal[(1-1):order - 1] + 1
        supdiag=diagonal[(2 + odd-1):2:order] - 1
        g=- ABCD[(supdiag-1)]
        c=np.array([- ABCD[0,(order + 3-1)],ABCD[(subdiag-1)]]).reshape(1,-1)
        a=np.zeros(1,order)
        for i in range(1,(order+1),2):
            a[(i-1)]=ABCD[(order + 1-1),(i-1)]
            ABCD[(order + 1-1),:]=ABCD[(order + 1-1),:] - a[(i-1)] * ABCD[(i-1),:]
        a[(2-1):2:order]=ABCD[(order + 1-1),(2-1):2:order]
        b=ABCD[(1-1):order + 1,(order + 2-1)].T
        for i in range(2,(order+1),2):
            b[(i-1)]=b[(i-1)] - c[(i-1)] * b[(i - 1-1)]
            if odd:
                b[(i-1)]=b[(i-1)] + g[(i / 2-1)] * b[(i + 1-1)]
        yscale=ABCD[(order + 2-1),(order + 1-1)]
        a=a * yscale
        b[(end-1)]=b[(end-1)] * yscale
    elif form in 'CIFF', 'Stratos':
        a = ABCD[(order + 1-1),(1-1):order]
        c=np.array([- ABCD[0,(order + 2-1)],ABCD[(subdiag[(1-1):end - 1]-1)]]).reshape(1,-1)
        g=- ABCD[(supdiag-1)]
        b=ABCD[:,(order + 1-1)].T
    else:
        fprintf(1,'%s error: Form %s is not yet supported.\\n',mfilename,form)
    return a,g,b,c
