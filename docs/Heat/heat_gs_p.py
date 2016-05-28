#coding: utf-8

L=0.1 #長さ
M=10 #分割数
dx=L/M #空間刻み
dt=1.0 #時間刻み
N=100 #時間ステップ数
alpha=80.2/(7874.0*440.0) #熱伝導率
gamma=dt/dx**2
a=alpha*gamma
T=[273.15 for i in range(M+1)]
EPS=1e-15 

def gaussseidel(T,preT):
	while True:
		d=0.0
		for i in range(1,M):
			n=(preT[i]+a*(T[i+1]+T[i-1]))/(1+2*a)
			d+=abs(n-T[i])
			T[i]=n
		if d < EPS:
			break

if __name__=='__main__':

	f = open('output','w')
	for i in range(1,N):
		T[0]=300.0
		T[M]=T[M-1]
		preT = list(T)
		gaussseidel(T,preT)
		f.write('# t=%ss\n'%(i*dt))
		for j in range(M+1):
			f.write('%s, %s\n'%(j*dx,T[j]))
		f.write('\n\n')
	f.close()
