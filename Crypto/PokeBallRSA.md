# Poke Ball RSA

## Challenge
We are given the value of n, e, c.
```
n = 498934084350094415783044823223130007435556803301613073259727203199325937230080661117917023582579699673759861892703348357714077684549303787581429366922208568924252052118455313229534699860304480039147103608782140303489222166267907007839021544433148286217133494762766492655602977085105487216032806292874190551319
e = 134901827939710543990222584187396847806193644190423846456160711527109836908087675183249532946675670587286594441908191054495871501233678465783530503352727362726294270065122447852357566161748618195216611965946646411519602447104878893524856862722902833460104389620397589021732407447981724307130484482495521398799
c = 100132888193232309251839777842498074992587507373917163874335385921940537055226546911990198769720313749286675018486390873216490470403470144298153410686092752282228631590006943913867497072931343354481759219425807850047083814816718302223434388744485547550941814186146959750515114700335721173624212499886218608818
```

As usual, I will start by running `rsactftool` (Reasoning: It is a RSA challenge).

After running `rsactftool` tool, it prints the following output:
```
e=2,c=9019127052844164572606928250741960583163943438936945828390420331200602392329
```

Sadly, the plaintext is not a flag. We need to continue to decrypt. Since the value of e is 2, it is a Rabin cryptosystem. Rabin cryptosystem, like RSA, uses asymmetric cryptography. 

Note: I will use the same `n` provided by the challenge (clarified with the author).

In short, Rabin's decryption (For more information: https://en.wikipedia.org/wiki/Rabin_cryptosystem):
1. Compute m<sub>q</sub> and m<sub>p</sub> (There are 2 different methods)
    1. Tonelli–shanks algorithm
    2. m<sub>p</sub> = c<sup>(p + 1) / 4</sup> mod p, m<sub>q</sub> = c<sup>(q + 1) / 4</sup> mod q (iff p = q = 3 (mod 4))
2. Extended euclidean algorithm
3. Chinese Remainder Theorem

In our case, either p or q is not 3 (mod 4), so I will perform Tonelli-shanks algorithm to compute m<sub>q</sub> and m<sub>p</sub>. Then, perform extended euclidean algorithm and lastly, perform chinese remainder theorem.

``` python
from Crypto.Util.number import *
import gmpy2
import itertools
from Crypto.Util.number import *
from functools import *

def ext_gcd(a: int, b: int):
    c0, c1 = a, b
    a0, a1 = 1, 0
    b0, b1 = 0, 1

    while c1 != 0:
        q, m = divmod(c0, c1)
        c0, c1 = c1, m
        a0, a1 = a1, (a0 - q * a1)
        b0, b1 = b1, (b0 - q * b1)

    return a0, b0, c0

n = 498934084350094415783044823223130007435556803301613073259727203199325937230080661117917023582579699673759861892703348357714077684549303787581429366922208568924252052118455313229534699860304480039147103608782140303489222166267907007839021544433148286217133494762766492655602977085105487216032806292874190551319
p = 16561429133925410152910558138426342813003607408793210510350820778117750488294158025334154782168295150504888994664227526968840978993964333891214657448288057
q = 30126269920030533619149872131139522468497774940789529153525410135192816439227647545595239038914424116121815731200170441720580775486313620726576435961039567

# obtained from rsactftool
e = 2
c = 9019127052844164572606928250741960583163943438936945828390420331200602392329

# tonelli–shanks algorithm
def legendre(a, p):
    return pow(a, (p - 1) // 2, p)
 
def tonelli(n, p):
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r

mp = tonelli(c,p)
mq = tonelli(c,q)

# Extended euclidean algorithm
yp, yq, _ = ext_gcd(p, q)

# Chinese Remainder Theorem
r = (yp * p * mq + yq * q * mp) % n
mr = n - int(r)
s = (yp * p * mq - yq * q * mp) % n
ms = n - int(s)
for num in [r,mr,s,ms]:
    print(long_to_bytes(num))
```

## Resource
https://medium.com/hackzone/hackzone-viii-ctf-rootcrypto-writeup-b1344f0d44d0
