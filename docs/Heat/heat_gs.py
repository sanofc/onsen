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
A=[[0.0 for i in range(M+1)] for j in range(M+1)]
EPS=1e-15 

def gaussseidel(T,preT,A):
	while True:
		d=0.0
		for i in range(M+1):
			sum=0.0
			for j in range(M+1):
				if i != j:
					sum+=(A[i][j]/A[i][i])*T[j]
			n=preT[i]/A[i][i]-sum
			d+=abs(n-T[i])
			T[i]=n
		if d < EPS:
			break

if __name__=='__main__':
	for i in range(1,M):
		A[i][i-1] = -a
		A[i][i] = 1+2*a
		A[i][i+1] = -a
	A[0][0]=1.0
	A[M][M]=1.0

	f = open('output','w')
	for i in range(1,N):
		T[0]=300.0
		T[M]=T[M-1]	
		preT = list(T)
		gaussseidel(T,preT,A)
		f.write('# t=%ss\n'%(i*dt))
		for j in range(M+1):
			f.write('%s, %s\n'%(j*dx,T[j]))
		f.write('\n\n')
	f.close()
