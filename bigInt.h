#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <math.h>

using namespace std;

const double PI = acos(-1);

class bigInt{
    private:
    vector<unsigned short>number;
    bool sign;

    bigInt(vector<unsigned short>num, bool sgn = true) {
        while(num.size() > 1 and num.back() == 0)
            num.pop_back();
        number = num;
        sign = sgn;
    }

    vector<complex<double>> fft(vector<complex<double>> a, bool invert = false) const {
        // from https://cp-algorithms.com/algebra/fft.html
        int n = a.size();
        if (n == 1)
            return vector<complex<double>>(1, a[0]);

        vector<complex<double>> a0(n / 2), a1(n / 2);
        for (int i = 0; 2 * i < n; i++) {
            a0[i] = a[2*i];
            a1[i] = a[2*i+1];
        }
        auto f0 = fft(a0, invert);
        auto f1 = fft(a1, invert);

        vector<complex<double>>f(a.size());

        double ang = 2 * PI / n * (invert ? -1 : 1);
        complex<double> w(1), wn(cos(ang), sin(ang));
        for (int i = 0; 2 * i < n; i++) 
        {
            f[i] = f0[i] + w * f1[i];
            f[i + n/2] = f0[i] - w * f1[i];
            if (invert) 
            {
                f[i] /= 2;
                f[i + n/2] /= 2;
            }   
            w *= wn;
        }
        return f;
    }

    
    public:

    bigInt(long long v = 0)
    {
        number = {};
        sign = (v > 0);
        if(v < 0)
        v = -v;
        if(v == 0)
        number = {0};
        while(v)
        {
            number.push_back(v % 10);
            v /= 10;
        }
    }

    bigInt(string s)
    {
        if(s[0]=='-')
            sign = 0;
        else
            sign = 1;
        number = {};
        for(int i=s.size()-1;i>=0;i--)
        {
            if(s[i] != '-')
            {
                number.push_back(s[i] - '0');
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const bigInt& obj) {
        if(!obj.sign)
        os << "-";
        for(int i = obj.number.size() - 1;i >= 0;i--)
        {
            os << (int)obj.number[i];
        }
        return os;
    }

    bigInt operator-() const {
        auto p = *this;
        p.sign = (!p.sign);
        return p;
    }

    bool operator<(const bigInt&other) const {
        if(sign)
        {
            if(!other.sign)
            {
                return false;
            }
            if(number.size() != other.number.size())
            {
                return number.size() < other.number.size();
            }
            for(int i = number.size() - 1;i >= 0;i--)
            {
                if(number[i] != other.number[i])
                {
                    return number[i] < other.number[i];
                }
            }
            return false;
        }
        else
        {
            if(other.sign)
            {
                return true;
            }
            return (-other) < (-*this);
        }
    }

    bool operator>(const bigInt&other) const {
        return (other) < (*this);
    }

    bool operator==(const bigInt&other) const {
        return !(*this < other) and !(other < *this);
    }

    bool operator>=(const bigInt&other) const {
        return !((*this) < (other));
    }

    bool operator<=(const bigInt&other) const {
        return !(other < (*this));
    }

    bigInt operator+(const bigInt& other) const {
        if(sign and other.sign)
        {
            unsigned short carry = 0;
            vector<unsigned short>result;
            for(int i = 0;i < max(number.size(), other.number.size());i++)
            {
                unsigned short a = 0, b = 0;
                if(i < number.size())
                    a = number[i];
                if(i < other.number.size())
                    b = other.number[i];
                unsigned short res = (a + b + carry);
                result.push_back(res % 10);
                carry = res/10;
            }   
            if(carry > 0)
                result.push_back(carry);
            return bigInt(result);
        }
        else
        if(sign and !other.sign)
        {
            return *this - (-other);
        }
        else
        if(!sign and other.sign)
        {
            return other - (-*this);
        }
        return - (-*this + (-other));
    }

    bigInt operator-(const bigInt& other) const {
        if(sign and other.sign)
        {
            if(*this >= other)
            {
                bool borrow = 0;
                vector<unsigned short>result;
                for(int i = 0;i < number.size();i++)
                {
                    int a = number[i];
                    int b = 0;
                    if(i < other.number.size())
                    {
                        b = other.number[i];
                    }
                    if(borrow)
                    {
                        if(a == 0)
                        {
                            a = 9;
                        }
                        else
                        {
                            borrow = false;
                            a--;
                        }
                    }
                    if(a < b)
                    {
                        a += 10;
                        borrow = true;
                    }
                    result.push_back(a - b);
                }
                return bigInt(result);
            }
            else
            {
                return -(other - *this);
            }
        }
        else
        if(sign and !other.sign)
        {
            return *this + (-other);
        }
        else
        if(!sign and other.sign)
        {
            return - (-*this + other);
        }
        return -other - (-*this);
    }

    bigInt operator*(const bigInt& other) const {

        bool sgn = (!(this->sign ^ other.sign));
        vector<complex<double>> a(this->number.begin(), this->number.end());
        vector<complex<double>> b(other.number.begin(), other.number.end());

        int n = 1;
        while(n < a.size() + b.size())
            n *= 2;
        a.resize(n);
        b.resize(n);

        auto fa = fft(a);
        auto fb = fft(b);

        for(int i = 0;i < n;i++)
        {
            fa[i]*=fb[i];
        }

        auto f = fft(fa, true);

        vector<unsigned short> ans;
        int carry = 0;
        for(int i = 0;i < n;i++)
        {
            int x = round(f[i].real());
            ans.push_back((carry + x)%10);
            carry = (carry + x)/10;
        }
        while(carry)
        {
            ans.push_back(carry % 10);
            carry /= 10;
        }
        return bigInt(ans, sgn);

    }

};

bigInt pow(bigInt a, int n)
{
    if(n==0)
    return 1;
    if(n%2==0)
    return pow(a*a, n/2);
    return a*pow(a, n-1);
}
