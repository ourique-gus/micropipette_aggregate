all:CalcForce

CalcForce:CalcForce.f95
	f2py3 --f90flags='-Wno-tab -fbounds-check -fcheck=all -fopenmp' -lgomp -c CalcForce.f95 -m CalcForce

clean:
	rm CalcForce*.so
