#include <iostream>
#include <cstdlib>

#define f(x,y) for (int x = 0; x < y; ++x)

using namespace std;

struct quat {
	long long q[4];
	
	quat() {
	}
	
	quat(long long a, long long b, long long c, long long d) {
		q[0] = a;
		q[1] = b;
		q[2] = c;
		q[3] = d;
	}
};

quat mult(quat x, quat y) {
	quat z;
	f(i,4) z.q[i] = x.q[0]*y.q[i]+x.q[i]*y.q[0];
	f(i,4) z.q[0] -= x.q[i]*y.q[i];
	f(i,3) z.q[i+1] += x.q[(i+1)%3+1]*y.q[(i+2)%3+1] - x.q[(i+2)%3+1]*y.q[(i+1)%3+1];
	f(i,4) z.q[i] /= 2;
	return z;
}

quat add(quat x, quat y) {
	quat z;
	f(i,4) z.q[i] = x.q[i] + y.q[i];
	return z;
}

quat sub(quat x, quat y) {
	quat z;
	f(i,4) z.q[i] = x.q[i] - y.q[i];
	return z;
}

void show(quat x) {
	if (x.q[0]%2 == 1) {
		f(i,3) cout << (x.q[i]/2.0) << " ";
		cout << (x.q[3]/2.0);
	} else {
		f(i,3) cout << (x.q[i]/2) << " ";
		cout << (x.q[3]/2);
	}
}

quat conj(quat x) {
	quat y;
	y.q[0] = x.q[0];
	f(i,3) y.q[i+1] = -x.q[i+1];
	return y;
}

long long norm(quat x) {
	long long nm = 0;
	f(i,4) nm += x.q[i]*x.q[i];
	return nm/4;
}

quat rquot(quat x, quat y) {
	quat ybx = mult(conj(y), x);
	long long ny = norm(y), rem[4], rnm = 0, r2;
	if (ny == 0) {
		cout << "Divided by zero!" << endl;
		f(i,4) ybx.q[i] = 0;
		return ybx;
	}
	f(i,4) {
		rem[i] = (ybx.q[i]%(2*ny)+2*ny)%(2*ny);
		if (rem[i] > ny) rem[i] -= 2*ny;
		r2 = rem[i] + ny;
		if (r2 > ny) r2 -= 2*ny;
		rnm += rem[i]*rem[i] - r2*r2;
	}
	if (rnm > 0) f(i,4) {
		if (rem[i] > 0) rem[i] -= ny;
		else rem[i] += ny;
	}
	f(i,4) ybx.q[i] = (ybx.q[i]-rem[i])/ny;
	return ybx;
}

quat gcld(quat x, quat y) {
	while (norm(y) != 0) {
		quat z = sub(x,mult(y,rquot(x,y)));
		x = y;
		y = z;
	}

	// post-processing: pick "nicest" right associate
	if (x.q[0]%2 != 0) {	// find a properly integral associate
		quat t = conj(x);
		f(i,4) {
			t.q[i] = (t.q[i]%4+4)%4;
			if (t.q[i] == 3) t.q[i] = -1;
		}
		x = mult(x,t);
	}
	long long nm = norm(x);
	if (nm%2 == 0) {	// just make it as positive as possible
		if (x.q[0] < 0) x = mult(x,quat(-2,0,0,0));
		f(i,3) if (x.q[i+1] < 0 && x.q[(i+1)%3+1] <= 0) {
			quat im(0,0,0,0);
			im.q[i+1] = 2;
			x = mult(x,im);
		}
		for (int i = 2; i >= 0; --i) if (x.q[i+1] < 0) {
			quat im(0,0,0,0);
			im.q[i+1] = 2;
			x = mult(x,im);
		}
	} else {	// choose an associate congruent to 1 mod 2
		bool one = (nm%4 == 1);
		f(i,4) if ((x.q[i]%4 == 0) ^ one) {
			quat im(0,0,0,0);
			im.q[i] = 2;
			x = mult(x,im);
			break;
		}
		if (x.q[0] < 0 || (x.q[0] == 0 && x.q[1] < 0)) x = mult(x,quat(-2,0,0,0));
	}
	return x;
}

int main(int argc, char **argv) {
	if (argc < 9) {
		cout << "Enter two (doubled) quaternions as arguments!" << endl;
		return 0;
	}
	quat a, b;
	bool good = true;
	f(i,4) {
		a.q[i] = atoi(argv[i+1]);
		if (i && (a.q[i]+a.q[0])%2 != 0) good = false;	// they didn't bother doubling the quaternion
	}
	if (!good) f(i,4) a.q[i] *= 2;
	good = true;
	f(i,4) {
		b.q[i] = atoi(argv[i+5]);
		if (i && (b.q[i]+b.q[0])%2 != 0) good = false;
	}
	if (!good) f(i,4) b.q[i] *= 2;
	
	cout << "(Product: ";
	show(mult(a,b));
	cout << ")\n" << endl;
	cout << "Greatest Common Left Divisor: ";
	quat l = gcld(a,b);
	show(l);
	cout << endl << "First Quotient: ";
	show(rquot(a,l));
	cout << endl << "Second Quotient: ";
	show(rquot(b,l));
	cout << endl;
	
	return 0;
}
