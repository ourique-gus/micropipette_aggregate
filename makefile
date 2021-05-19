all:CalcForce

CalcForce:CalcForce.f95
	f2py --f90flags='-Wno-tab -fbounds-check -fcheck=all' -lgomp -c CalcForce.f95 -m CalcForce
	f2py --f90flags='-Wno-tab -fbounds-check -fcheck=all -fopenmp' -lgomp -c CalcForce.f95 -m CalcForceOMP

clean:
	rm CalcForce*.so
