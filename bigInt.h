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
    vector<int>number;
    bool sign;
    int base = 2;

    bigInt(vector<int>num, bool sgn = true, int base = 10) {
        while(num.size() > 1 and num.back() == 0)
            num.pop_back();
        number = num;
        sign = sgn;
        if(*this == 0)
            sign = 1;
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
        // TODO: negatives
        if(other == 0)
        {
            return {1 / 0, 1 / 0};
        }

        bigInt quotient = 0, remainder = 0;
        for(int i = number.size() - 1;i >= 0;i--)
        {
            remainder = (remainder << 1);
            quotient = (quotient << 1);
            remainder = remainder + number[i];
            
            if(remainder >= other)
            {
                remainder = remainder - other;
                quotient = quotient + 1;
            }
        }
        return {quotient, remainder};
    }

    
    public:

    bigInt(long long v = 0)
    {
        number = {};
        sign = (v >= 0);
        if(v < 0)
        v = -v;
        if(v == 0)
        number = {0};
        while(v)
        {
            number.push_back(v % base);
            v /= base;
        }
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
        vector<int>zero = {0};
        if(p.number == zero)
        return p;
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
        if(sign and other.sign)
        {
            int carry = 0;
            vector<int>result;
            for(int i = 0;i < max(number.size(), other.number.size());i++)
            {
                int a = 0, b = 0;
                if(i < number.size())
                    a = number[i];
                if(i < other.number.size())
                    b = other.number[i];
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
                vector<int>result;
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

    bigInt operator&(const bigInt& other) const {
        vector<int>res;
        for(int i = 0;i < min(other.number.size(), number.size());i++)
        {
            res.push_back(number[i] & other.number[i]);            
        }
        return bigInt(res);
    }
    
    bigInt operator|(const bigInt& other) const {
        vector<int>res;
        for(int i = 0;i < max(other.number.size(), number.size());i++)
        {
            int x = 0, y = 0;
            if(i < number.size())
            x = number[i]; 
            if(i < other.number.size())
            y = other.number[i];   
            res.push_back(x | y); 
        }
        return bigInt(res);
    }

    bigInt operator^(const bigInt& other) const {
        vector<int>res;
        for(int i = 0;i < max(other.number.size(), number.size());i++)
        {
            int x = 0, y = 0;
            if(i < number.size())
            x = number[i]; 
            if(i < other.number.size())
            y = other.number[i];   
            res.push_back(x ^ y); 
        }
        return bigInt(res);
    }

    bigInt operator>>(int p) const {
        auto res = number;
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
        auto res = number;
        reverse(res.begin(), res.end());
        for(int i = 0;i < p;i++)
        {
            res.push_back(0);
        }
        reverse(res.begin(), res.end());
        return bigInt(res);
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

    }

    bigInt operator%(const bigInt& other) const {

        return divide(other).second;

    }

    int number_of_bits() const{
        return number.size();
    }

};

bigInt pow(bigInt a, int n)
{
    if(n == 0)
        return 1;
    if(n % 2 == 0)
        return pow(a * a, n / 2);
    return a * pow(a, n - 1);

}

bigInt gcd(bigInt a, bigInt b)
{
    if(b == 0)
    return a;
    return gcd(b, a % b);
}

#endif 
