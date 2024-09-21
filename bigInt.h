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
    int bits = 1;
    long long base = (1 << bits);

    bigInt(vector<int>num, bool sgn = false) {
        while(num.size() > 1 and num.back() == 0)
            num.pop_back();
        digits = num;
        sign = sgn;
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
        for(int i = 0;i < p / bits and res.size();i++)
        {
            res.pop_back();
        }
        if(res.size() == 0)
        return bigInt(0);
        
        int k = p % bits;
        res[res.size() - 1] >>= k;

        for(int i = res.size() - 2;i >= 0;i--)
        {
            res[i + 1] |= ((res[i] & ((1<<k) - 1))<<(bits - k));
            res[i] >>= k;
        }
        reverse(res.begin(), res.end());
        return bigInt(res);
    }

    bigInt operator<<(int p) const {
        auto res = digits;
        res.push_back(0);
        reverse(res.begin(), res.end());
        int k = p % bits;
        for(int i = 1;i < res.size();i++)
        {
            res[i - 1] |= (res[i] >> (bits - k));
            res[i] <<= k;
            res[i] = (res[i] & ((1 << bits) - 1));
        }
        for(int i = 0;i < p / bits;i++)
        {
            res.push_back(0);
        }
        reverse(res.begin(), res.end());
        return bigInt(res);
    }

    bigInt operator*(const bigInt& other) const {

        bool sgn = ((sign ^ other.sign));
        vector<complex<double>> a(digits.begin(), digits.end());
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
        long long carry = 0;
        for(int i = 0;i < n;i++)
        {
            long long x = round(f[i].real());
            ans.push_back(1ll * (carry + x) % base);
            carry = 1ll * (carry + x) / base;
        }

        while(carry)
        {
            ans.push_back(1ll * carry % base);
            carry /= base;
        }
        return bigInt(ans, sgn);

    }

    bigInt operator/(const bigInt& other) const {

        bigInt div = *this;

        if(div < other)
        {
            return 0;
        }


        int n = div.number_of_bits() + other.number_of_bits();
        bigInt power_of_two = (bigInt(1) << (n + 1));

        bigInt x = div - other;
        bigInt last = -1, before_last = -1;// two last values of x

        while(x != last and x != before_last)
        {
            cerr << x.number_of_bits()<<endl;
            before_last = last;
            last = x;
            x = ((x * (power_of_two - x * other)) >> n);
        }
        bigInt ans = ((div * x) >> n);

        if(div - (ans * other) >= other)
        {
            ans = (ans + bigInt(1));
        }

        return ans;

    }

    bigInt operator%(const bigInt& other) const {
        bigInt div = *this;
        auto x = (div / other);
        return div - other * x;
    }

    vector<int> change_base(int _base = 10) const {

        bigInt x = *this;

        if(x == bigInt(0))
        {
            return {0};
        }

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
        int x = digits.back();
        int c = 0;
        while(x)
        {
            c++;
            x >>= 1;
        }
        return bits * (digits.size() - 1) + c;
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
