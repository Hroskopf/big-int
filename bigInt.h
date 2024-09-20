#ifndef BIGINT
#define BIGINT

#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <math.h>

using namespace std;

class bigInt{
    private: 
    vector<int>digits;
    bool sign;
    int base = 2;

    bigInt(vector<int>num, bool sgn = false, int _base = 2) {
        while(num.size() > 1 and num.back() == 0)
            num.pop_back();
        digits = num;
        sign = sgn;
        base = _base;
        if(*this == (bigInt)0)
            sign = 0;
    }

    vector<complex<double>> fft(vector<complex<double>> a, bool invert = false) const {
        // inspired by https://cp-algorithms.com/algebra/fft.html
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

        const double PI = acos(-1);

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

    pair<bigInt, bigInt> divide(const bigInt& other) const {
        // O(n ^ 2) -- slow
        // TODO: faster
        if(other == (bigInt)0)
        {
            return {1 / 0, 0};
        }

        bigInt quotient = 0, remainder = 0;
        for(int i = digits.size() - 1;i >= 0;i--)
        {
            remainder = (remainder << 1);
            quotient = (quotient << 1);
            remainder = remainder + (bigInt)digits[i];
            
            if(remainder >= other)
            {
                remainder = remainder - other;
                quotient = quotient + (bigInt)1;
            }
        }
        return {quotient, remainder};
    }

    
    public:

    bigInt(long long v = 0)
    {
        digits = {};
        sign = (v < 0);
        if(v < 0)
        v = -v;
        if(v == 0)
        digits = {0};
        while(v)
        {
            digits.push_back(v % base);
            v /= base;
        }
    }

    operator int() const {

        int res = 0;
        int power_of_base = 1;
        for(int i = 0;i < digits.size();i++)
        {
            res += digits[i] * power_of_base;
            power_of_base = power_of_base * base;
        }
        if(sign)
        res = -res;
        return res;
    }

    bigInt(string s)
    {
        // TODO
        // if(s[0]=='-')
        //     sign = 0;
        // else
        //     sign = 1;
        // number = {};
        // base = 10;
        // for(int i=s.size()-1;i>=0;i--)
        // {
        //     if(s[i] != '-')
        //     {
        //         number.push_back(s[i] - '0');
        //     }
        // }
    }

    friend std::ostream& operator<<(std::ostream& os, const bigInt& obj) {
        string s = obj.to_decimal();
        os << s;
        return os;
    }

    bigInt operator-() const {
        auto p = *this;
        vector<int>zero = {0};
        if(p.digits == zero)
        return p;
        p.sign = (!p.sign);
        return p;
    }

    bool operator<(const bigInt&other) const {
        if(!sign)
        {
            if(other.sign)
            {
                return false;
            }
            if(digits.size() != other.digits.size())
            {
                return digits.size() < other.digits.size();
            }
            for(int i = digits.size() - 1;i >= 0;i--)
            {
                if(digits[i] != other.digits[i])
                {
                    return digits[i] < other.digits[i];
                }
            }
            return false;
        }
        else
        {
            if(!other.sign)
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

    bool operator!=(const bigInt&other) const {
        return !(*this == other);
    }

    bool operator>=(const bigInt&other) const {
        return !((*this) < (other));
    }

    bool operator<=(const bigInt&other) const {
        return !(other < (*this));
    }

    bigInt operator+(const bigInt& other) const {
        if(!sign and !other.sign)
        {
            int carry = 0;
            vector<int>result;
            for(int i = 0;i < max(digits.size(), other.digits.size());i++)
            {
                int a = 0, b = 0;
                if(i < digits.size())
                    a = digits[i];
                if(i < other.digits.size())
                    b = other.digits[i];
                int res = (a + b + carry);
                result.push_back(res % base);
                carry = res / base;
            }   
            while(carry > 0)
            {
                result.push_back(carry % base);
                carry /= base;
            }

            return bigInt(result);
        }
        else
        if(!sign and other.sign)
        {
            return *this - (-other);
        }
        else
        if(sign and !other.sign)
        {
            return other - (-*this);
        }
        return - (-*this + (-other));
    }

    bigInt operator-(const bigInt& other) const {
        if(!sign and !other.sign)
        {
            if(*this >= other)
            {
                bool borrow = 0;
                vector<int>result;
                for(int i = 0;i < digits.size();i++)
                {
                    int a = digits[i];
                    int b = 0;
                    if(i < other.digits.size())
                    {
                        b = other.digits[i];
                    }
                    if(borrow)
                    {
                        if(a == 0)
                        {
                            a = base - 1;
                        }
                        else
                        {
                            borrow = false;
                            a--;
                        }
                    }
                    if(a < b)
                    {
                        a += base;
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
        if(!sign and other.sign)
        {
            return *this + (-other);
        }
        else
        if(sign and !other.sign)
        {
            return - (-*this + other);
        }
        return -other - (-*this);
    }

    bigInt operator&(const bigInt& other) const {
        vector<int>res;
        for(int i = 0;i < min(other.digits.size(), digits.size());i++)
        {
            res.push_back(digits[i] & other.digits[i]);            
        }
        return bigInt(res);
    }
    
    bigInt operator|(const bigInt& other) const {
        vector<int>res;
        for(int i = 0;i < max(other.digits.size(), digits.size());i++)
        {
            int x = 0, y = 0;
            if(i < digits.size())
            x = digits[i]; 
            if(i < other.digits.size())
            y = other.digits[i];   
            res.push_back(x | y); 
        }
        return bigInt(res);
    }

    bigInt operator^(const bigInt& other) const {
        vector<int>res;
        for(int i = 0;i < max(other.digits.size(), digits.size());i++)
        {
            int x = 0, y = 0;
            if(i < digits.size())
            x = digits[i]; 
            if(i < other.digits.size())
            y = other.digits[i];   
            res.push_back(x ^ y); 
        }
        return bigInt(res);
    }

    bigInt operator>>(int p) const {
        auto res = digits;
        reverse(res.begin(), res.end());
        for(int i = 0;i < p and res.size();i++)
        {
            res.pop_back();
        }
        if(res.size() == 0)
        res = {0};
        reverse(res.begin(), res.end());
        return bigInt(res);
    }

    bigInt operator<<(int p) const {
        auto res = digits;
        reverse(res.begin(), res.end());
        for(int i = 0;i < p;i++)
        {
            res.push_back(0);
        }
        reverse(res.begin(), res.end());
        return bigInt(res);
    }

    bigInt operator*(const bigInt& other) const {


        bool sgn = ((this->sign ^ other.sign));
        vector<complex<double>> a(this->digits.begin(), this->digits.end());
        vector<complex<double>> b(other.digits.begin(), other.digits.end());

        int n = 1;
        while(n < a.size() + b.size())
            n *= 2;
        a.resize(n);
        b.resize(n);

        auto fa = fft(a);
        auto fb = fft(b);

        for(int i = 0;i < n;i++)
        {
            fa[i] *= fb[i];
        }

        auto f = fft(fa, true);

        vector <int> ans;
        int carry = 0;
        for(int i = 0;i < n;i++)
        {
            int x = round(f[i].real());
            ans.push_back((carry + x) % base);
            carry = (carry + x) / base;
        }
        while(carry)
        {
            ans.push_back(carry % base);
            carry /= base;
        }
        return bigInt(ans, sgn);

    }

    bigInt operator/(const bigInt& other) const {

        return divide(other).first;

        // if(*this < other)
        // {
        //     return 0;
        // }

        // int n = digits.size() + 5;
        

    }

    bigInt operator%(const bigInt& other) const {

        return divide(other).second;

    }

    vector<int> change_base(int _base = 10) const {

        bigInt x = *this;
        bool sgn = x.sign;
        vector<int>res;
        bigInt b = _base;
        while(x != (bigInt)0)
        {
            res.push_back((x % b));
            x = x / b;
        }
        return res;

    }

    int number_of_bits() const{
        return digits.size();
    }

    string to_bin() const {
        string s = "";
        if(sign)s += '-';
        for(int i = digits.size() - 1;i>=0;i--)
        {
            s += ('0' + digits[i]);
        }
        return s;

    }

    string to_hex() const {
        auto d = change_base(16);
        string s = "";
        if(sign)s += '-';
        for(int i = d.size() - 1;i>=0;i--)
        {
            if(d[i] < 10)
            s += ('0' + d[i]);
            else
            s += ('A' + d[i] - 10);
        }
        return s;

    }

    string to_decimal () const {
        auto d = change_base(10);
        string s = "";
        if(sign)s += '-';
        for(int i = d.size() - 1;i>=0;i--)
        {
            if(d[i] < 10)
            s += ('0' + d[i]);
        }
        return s;

    }

};

bigInt pow(bigInt a, int n)
{
    if(n == 0)
        return 1;
    if(n % 2 == 0)
        return pow(a * a, n / 2);
    return a * pow(a * a, n / 2);

}

bigInt gcd(bigInt a, bigInt b)
{
    if(b == (bigInt)0)
    return a;
    return gcd(b, a % b);
}

#endif 
