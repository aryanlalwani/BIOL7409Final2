:Reference :Colbert and Pan 2002

NEURON	{
	SUFFIX nav11
	USEION na READ ena WRITE ina
	RANGE gbar, g, ina, vtha, vthi, vthis, qa, qi, qs
	RANGE m, n, s, stau_scale, htau_scale 
	RANGE sTau, hTau
	RANGE vShift
	GLOBAL q10
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	vShift = 0 (mV) : try -8 mV to match Schmidt-Hieber 2010
	gbar = 0.2 (S/cm2)
	vtha = -21.2 (mV) : v 1/2 for activation (m) - Table 1 Spampanato 2004 Nav1.1 channel
	vthi = -39.7 (mV) : v 1/2 for inactivation (h)
    :vthis = -46.1 (mV) : v 1/2 for slow inactivation (s) - 2004 paper 
	vthis = -46.0 (mV) : v 1/2 for slow inactivation (s) - 2005 paper 
	qa = 4.9 (1) : activation slope
	qi = 7.7 (1) : inactivation slope	    
    	qs = 5.4 (1) : slow inactivation slow - 2004 paper
     :qs = 6.6 (1) : slow inactivation slow - 2005 paper
	q10  = 2.3			: temperature sensitivity
    stau_scale = 1 (1) : scaling factor for slow inactivation (sTau)
    htau_scale = 1 (1) : scaling factor for fast inactivation (hTau)
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	g	(S/cm2)
	mInf (1)
	mTau (ms)
	hInf (1)
	hTau (ms)	
	sInf (1)
	sTau (ms)
	celsius		(degC)
}

STATE	{
	m
	h
    s
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g = gbar*m*m*m*h*s
	ina = g*(v-ena)
}

DERIVATIVE states	{
	rates(q10)
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
	s' = (sInf-s)/sTau
}

INITIAL{
	rates(q10)
	m = mInf
	h = hInf
	s = sInf
}

PROCEDURE rates(q10){
  LOCAL qt
  :qt = 2.3^((34-21)/10)
  qt = q10^((celsius-21)/10)	
  UNITSOFF
    if(v-vShift == vtha){
    	v = v+0.0001
    }		
	mInf = 1/(1 + exp(-((v-vShift) - vtha)/qa))
	mTau = 0.15/qt
	
    if(v-vShift == vthi){
      v = v + 0.0001
    }		

	hInf = 1/(1 + exp(((v-vShift) - vthi)/qi))
	hTau = htau_scale*20.1*exp(-0.5*(((v-vShift) - -61.4)/32.7)^2)/qt		

	if (v-vShift == vthis){
      v = v + 0.0001
    }
	:sInf = 1/(1 + exp(((v-vShift) - vthis)/qs))
	:sTau = stau_scale*1000*(106.7*exp(-0.5*(((v-vShift) - -52.7)/18.3)^2))/qt		
	sInf = 1/(1+exp(((v-vShift) - vthis)/qs))		
	sTau = stau_scale*1000*(140.4*exp(-0.5*((v+71.3)/30.9)^2))/qt
	
	UNITSON
}
