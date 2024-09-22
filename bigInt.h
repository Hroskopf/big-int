#ifndef BIGINT
#define BIGINT

#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <math.h>

using namespace std;

class bigInt;

bigInt pow(bigInt a, int n);

bigInt convert(vector<int>num, int _base);


class bigInt{
    public:
    vector<int>digits;
    bool sign;
    

    private: 

    int bits = 10;
    int base = (1 << bits);
    vector<complex<double>> fft(vector<complex<double>> a, bool invert = false) const {

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

    vector<int>change_base(int _base, int order) const{

        bigInt x = *this;
        bigInt b = _base;

        if(x < bigInt(0))
        x = -x;
        
        if(order == 0)
        {
            return x.digits;
        }
        
        bigInt c = pow(b, order);

        bigInt l = x % c, r = x / c;
        vector<int>ans = l.change_base(_base, order / 2);
        for(int i:r.change_base(_base, order / 2))
            ans.push_back(i);
        return ans;
    }

    bigInt(vector<int>num, bool sgn = false) {

        while(num.size() > 1 and num.back() == 0)
            num.pop_back();
        digits = num;
        sign = sgn;
        if(*this == (bigInt)0)
            sign = 0;
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
        vector<int>num;
        int f = 0;
        if(s[0] == '-')
        {
            sign = 1;
            f++;
        }
        for(int i = f;i < s.size();i++)num.push_back(s[i] - '0');
        bigInt x = convert(num, 10);
        digits = x.digits;
        sign = (s[0] == '-');
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
        // for(auto i:f)cerr<<i;cerr<<endl;
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

    vector<int> change_base(int _base = 10) const {
        
        int order = 1;
        bigInt b = _base;        
        bigInt x = *this;
        if(x < bigInt(0))
        x = -x;
        while(pow(b, order + 1) <= x)
        {
            order = order * 2;
        }

        auto res = change_base(_base, order);
        while(res.size() > 1 and res.back() == 0)res.pop_back();
        return res;

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
        if(sign)s = '-';
        for(int i = d.size() - 1;i>=0;i--)
        {
            if(d[i] < 10)
                s += ('0' + d[i]);
        }
        return s;

    }

};

inline bigInt convert_rec(vector<int>num, bigInt _base, int order)
{
    if(order == 1)
    {
        return bigInt(num[0]);
    }
    vector<int>l, r;
    for(int i = 0;i < order / 2;i++)
    l.push_back(num[i]);
    for(int i = order / 2;i < order;i++)
    r.push_back(num[i]);
    bigInt L = convert_rec(l, _base, order / 2);
    bigInt R = convert_rec(r, _base, order / 2);
    return R + L * pow(_base, order / 2);
}

bigInt convert(vector<int>num, int _base = 10)
{
    // num is vector of digits of the number in given base. Most most significant bite first
    int order = 1;
    while(num.size() > order)
    {
        order *= 2;
    }
    reverse(num.begin(), num.end());
    num.resize(order);
    reverse(num.begin(), num.end());
    return convert_rec(num, bigInt(_base), order);
}

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
