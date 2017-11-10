iterations = 1000

pop = [6000 1000 1000 0 0 0 0
2000 5000 1000 0 10 10 0
500 3000 12000 0 30 0 0
0 0 10 0 0 2 0
2 10 20 10 3 20 0
5 2 0 0 0 10000 2000
3 1 0 0 1000 1000 5000].*1000

road = [11 122122
12 11132122
13 1223
21 1112223331
22 111221233233
23 13223233
31 2132
32 31223323
33 3221321223
35 25
25 35
43 5354
53 4354
52 61
61 5271
75 76
76 666777
66 766777
67 667677
77 666776]

aeroplane = [11 33
33 1177
77 33]

port = [33 43
33 35
43 33
35 33
53 52
52 53
25 26
26 25
26 46
46 26
34 56
56 34
56 66
66 56]

L = vec([.6 .5 .5 0 0 0 0
.5 .6 .6 0 .4 .3 0
.5 .6 .7 0 .3 0 0
0 0 .5 0 0 .3 0
.4 .5 .3 .2 .4 .5 0
.5 .3 0 0 0 .9 .8
.6 .4 0 0 .8 .8 .9]')'

beta = 0.7

B = vec([10000 1000 1000 0 0 0 0
1500 8000 1000 0 200 200 0
500 2000 20000 0 200 0 0
0 0 200 0 0 20 0
40 200 400 200 30 300 0
80 40 0 0 0 10000 1200
50 20 0 0 600 900 6000]')'

pop_a = vec(pop')'
pop_all = pop_a
for uk in range(1,size(pop)[1]^2-1)
  pop_all = vcat(pop_all,pop_a)
end

for jkl in range(1,2*size(pop)[1]^2)
  B = hcat(B,0)
end

Lst = copy(L)
for jkl in range(1,2*size(pop)[1]^2)
  Lst = hcat(Lst,0)
end

D = B.*(1.-Lst)

#----------------------------------------#
# Roads

roadmax = size(road)[1]

roads = eye(size(pop)[1]^2)
ports = eye(size(pop)[1]^2)
aeroplanes = eye(size(pop)[1]^2)

roadstart = 1
while roadstart <= roadmax
  place = road[roadstart,1]
  pdg = digits(place)
  pdg2 = pdg[end:-1:1,end:-1:1]
  down = (pdg2[1]-1)*7+pdg2[2]

  other = road[roadstart,2]
  digs = digits(other)[end:-1:1,end:-1:1]
  len = size(digs)[1]
  for xyz in range(1,convert(Int64,len/2))
    down2 = (digs[(xyz*2)-1]-1)*size(pop)[1]+digs[(xyz*2)]
    roads[down,down2]=1
  end
  roadstart += 1
end
currentroad = (1-eye(size(pop)[1]^2)).*roads.*pop_all
totalroad = sum(currentroad,2)
roadx = 1.6 /(1.6+0.4+0.1)
ratioroad = currentroad ./ totalroad

roadcorrect = copy(ratioroad)
for m=1:size(roadcorrect,1)
	for l=1:size(roadcorrect,2)
  		 isnan(roadcorrect[m,l]) ? roadcorrect[m,l]=0 :	 roadcorrect[m,l]=roadcorrect[m,l]
	end
end

#----------------------------------------#
# Aeroplane

aeromax = size(aeroplane)[1]
aerostart = 1
while aerostart <= aeromax
  place = aeroplane[aerostart,1]
  pdg = digits(place)
  pdg2 = pdg[end:-1:1,end:-1:1]
  down = (pdg2[1]-1)*7+pdg2[2]

  other = aeroplane[aerostart,2]
  digs = digits(other)[end:-1:1,end:-1:1]
  len = size(digs)[1]
  for xyz in range(1,convert(Int64,len/2))
    down2 = (digs[(xyz*2)-1]-1)*size(pop)[1]+digs[(xyz*2)]
    aeroplanes[down,down2]=1
  end
  aerostart += 1
end
currentaero = (1-eye(size(pop)[1]^2)).*aeroplanes.*pop_all
totalaero = sum(currentaero,2)
aerox = 0.1 /(1.6+0.4+0.1)
ratioaero = currentaero ./ totalaero

aerocorrect = copy(ratioaero)
for m=1:size(aerocorrect,1)
	for l=1:size(aerocorrect,2)
  		 isnan(aerocorrect[m,l]) ? aerocorrect[m,l]=0 :	 aerocorrect[m,l]=aerocorrect[m,l]
	end
end

#----------------------------------------#
# Ports

portmax = size(port)[1]

portstart = 1
while portstart <= portmax
  place = port[portstart,1]
  pdg = digits(place)
  pdg2 = pdg[end:-1:1,end:-1:1]
  down = (pdg2[1]-1)*7+pdg2[2]

  other = port[portstart,2]
  digs = digits(other)[end:-1:1,end:-1:1]
  len = size(digs)[1]
  for xyz in range(1,convert(Int64,len/2))
    down2 = (digs[(xyz*2)-1]-1)*size(pop)[1]+digs[(xyz*2)]
    ports[down,down2]=1
  end
  portstart += 1
end
currentport = (1-eye(size(pop)[1]^2)).*ports.*pop_all
totalport = sum(currentport,2)
portx = 0.4 /(1.6+0.4+0.1)
ratioport = currentport ./ totalport

portcorrect = copy(ratioport)
for m=1:size(portcorrect,1)
	for l=1:size(portcorrect,2)
  		 isnan(portcorrect[m,l]) ? portcorrect[m,l]=0 :	 portcorrect[m,l]=portcorrect[m,l]
	end
end

#----------------------------------------#
## ADDITIONS
rr = roadcorrect.*roadx.*0.015
pp = portcorrect.*portx.*0.015
aa = aerocorrect.*aerox.*0.015
needleft = 1 .- sum(rr+pp+aa,2)

eyeget = eye(size(pop)[1]^2) .* needleft

ssss = eyeget + aa + rr + pp

#----------------------------------------#

large = size(ssss)[1]

Ipart = zeros(1,size(pop)[1]^2)
Ipart[1,13]=1
Rpart = zeros(1,size(pop)[1]^2)

SIR = hcat(pop_a,Ipart,Rpart)

table = ["s1"]
for x in range(2,size(pop)[1]^2-1)
  h = string("s",x)
  table = hcat(table,h)
end
for x in range(1,size(pop)[1]^2)
  h = string("i",x)
  table = hcat(table,h)
end
for x in range(1,size(pop)[1]^2)
  h = string("r",x)
  table = hcat(table,h)
end

old_SIR = zeros(1,large*3)

mosquito = zeros(1,large)
old_ir = zeros(large,1)

#----------------------------------------#

for x in range(1,iterations)

  temp = rand(10:44,1,large)
  prec = rand(100:400,1,large)

#----------------------------------------#
  # Temperature first
  sd = 7;  mean = 27;  sc = 13.5
  h = (1/(sd*sqrt(2*pi)))*sc
  f(t) = (exp((-(t-mean).^2)/(2*sd^2)))*h

  T=f(temp)
  total = sum(T)
  ss = size(T)[2]
  temperature = eye(ss)

  for uu in 1:ss
    Q = copy(vec(T))
    Q[uu] = 0
    jj = 0

    for x in Q
      if x == 0
        d = T[uu]
      else
        y = total-T[uu]
        b = x/y
        c = 1-T[uu]
        d = c*b
      end
      jj += 1

      temperature[uu,jj] = d
    end
  end
#----------------------------------------#

  # Precipitation next
  sd = 50;  mean = 250;  sc = 100
  h = (1/(sd*sqrt(2*pi)))*sc
  f(t) = (exp((-(t-mean).^2)/(2*sd^2)))*h

  T=f(prec)
  total = sum(T)
  ss = size(T)[2]
  precipitation = eye(ss)

  for uu in 1:ss
    Q = copy(vec(T))
    Q[uu] = 0
    jj = 0

    for x in Q
      if x == 0
        d = T[uu]
      else
        y = total-T[uu]
        b = x/y
        c = 1-T[uu]
        d = c*b
      end
      jj += 1

      precipitation[uu,jj] = d
    end
  end
#----------------------------------------#

  # Living standards next
  T = 1 - L
  total = sum(T)
  ss = size(T)[2]
  living = eye(ss)

  for uu in 1:ss
    Q = copy(vec(T))
    Q[uu] = 0
    jj = 0

    for x in Q
      if x == 0
        d = T[uu]
      else
        y = total-T[uu]
        b = x/y
        c = 1-T[uu]
        d = c*b
      end
      jj += 1

      living[uu,jj] = d
    end
  end
#----------------------------------------#

  # Now probability for mosquito

  mosq_p = temperature*precipitation*living

  mosq_b = rand(1000:2000,1,large)
  mosq_d = rand(500:600,1,large)
  w_scalar = 0.7

  mosquito = (mosquito + (mosq_b*0.7 - mosq_d))*mosq_p

  # Now scaling for mosquito

  addup = sum(mosquito)
  elements = size(mosquito)[2]
  average = addup / elements

  addon_all = 0
  for ippa in mosquito
    addon = (ippa - average)^2
    addon_all = addon + addon_all
  end

  sd = (1/elements*addon_all)^0.5
  z = (mosquito .- average) ./ sd

  a = 111.*atan(27.*z/294)
  b = -358.*z./23
  c = b+a
  d = exp(c)
  attack_rate = 1./(1.+d)
#----------------------------------------#

  # S -> I matrix

  u = SIR[large+1:2*large] ./ (SIR[1:large]+SIR[large+1:2*large]+SIR[2*large+1:3*large])
  # pq = (SIR[large+1:2*large] .- old_SIR[1+large:2*large])
  # qrd = (SIR[1:large] .- old_SIR[1:large])

  ratio = 0.5
  IR = attack_rate' .* u .* ratio
  # constant = 0.1
  # IRRR = (2*exp(1/2*constant*x))/(1+exp(1/2*constant*x))
  # IRR = IRRR .* (SIR[10:18]./(SIR[1:9].+SIR[10:18].+SIR[19:27]))
  # IR = IRRR .* IRR
  kk = 1
  infection = zeros(elements)
  for xp in isnan(IR)*1
    if xp == 0
      infection[kk] = IR[kk]
    else
      infection[kk] = 0.0
    end
    kk += 1
  end

  old_SIR = copy(SIR)
  old_ir = copy(IR)

  s_i = infection' .* eye(elements)
  global(s_i)

#----------------------------------------#

  # I -> R matrix
  opp = 20
  i_r1 = L ./opp
  i_r = eye(elements).*i_r1

#----------------------------------------#

  # S -> S matrix
  nrow = size(ssss)[1]
  h = sum(eye(nrow,nrow)-s_i,2)

  s_s = ssss.*h

#----------------------------------------#

  # I -> I matrix
  leftover = (1-sum(i_r,2))
  scaled = leftover.*s_s
  other = (scaled.*(1-eye(large))).*0.5
  othertotal = sum(other,2)
  correct = (scaled.+othertotal).*eye(large).+other
  i_i = (leftover-sum(correct,2).+correct).*eye(large).+correct.*(1-eye(large))


#----------------------------------------#
  # R -> R matrix

  gg = (leftover.-s_s).*eye(nrow,nrow).*2/nrow
  hh = 1 .- sum(s_s+gg,2)
  qp = (1.-eye(nrow,nrow)).*s_s
  vv = qp ./ sum(qp,2) .*hh

  ooo = 1
  for xpd in isnan(vv)*1
    if xpd == 1
      vv[ooo] = 0
    end
    ooo += 1
  end
  r_r = gg+s_s+vv

#----------------------------------------#

  z = zeros(nrow,nrow)
  P1 = hcat(s_s,s_i,z)
  P2 = hcat(z,i_i,i_r)
  P3 = hcat(z,z,r_r)

  p = vcat(P1,P2,P3)

  SIR = (SIR+(B*beta-D))*p
  table = vcat(table,SIR)

end

#----------------------------------------#

table = table'
starting = ["Date"]
for x in range(1,iterations)
  starting = hcat(starting,x)
end

table = vcat(starting,table)
writecsv("Data.csv",table, header = false)
