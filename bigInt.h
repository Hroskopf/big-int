#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class bigInt{
    private:
    vector<unsigned short>number;
    bool sign;

    bigInt(vector<unsigned short>num, bool sgn = true)
    {
        while(num.size() > 1 and num.back() == 0)
            num.pop_back();
        number = num;
        sign = sgn;
    }

    bigInt multiple_by_ten_power(bigInt a, int n) const
    {
        if(a == 0)
        return a;
        auto b = a;
        reverse(b.number.begin(), b.number.end());
        for(int i=0;i<n;i++)
        {
            b.number.push_back(0);
        }
        reverse(b.number.begin(), b.number.end());
        return b;
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
        auto a = *this;
        auto b = other;
        if(max(a.number.size(), b.number.size()) == 1)
        {
            int x = a.number[0];
            int y = b.number[0];
            if(!a.sign)
            x = -x;
            if(!b.sign)
            y = -y;
            return bigInt(x*y);
        }
        vector<unsigned short>A, B, C, D;
        int n = max(a.number.size(), b.number.size()) / 2;
        for(int i=0;i<max(a.number.size(), b.number.size());i++)
        {
            if(i<n)
            {
                if(i < a.number.size())
                B.push_back(a.number[i]);
                else
                B.push_back(0);
                if(i < b.number.size())
                D.push_back(b.number[i]);
                else
                D.push_back(0);
            }
            else
            {
                if(i < a.number.size())
                A.push_back(a.number[i]);
                else
                A.push_back(0);
                if(i < b.number.size())
                C.push_back(b.number[i]);
                else
                C.push_back(0);
            }
        }
        bigInt A1 = bigInt(A, a.sign);
        bigInt B1 = bigInt(B, a.sign);
        bigInt C1 = bigInt(C, b.sign);
        bigInt D1 = bigInt(D, b.sign);
        bigInt X = (A1 + B1) * (C1 + D1);
        bigInt Y = A1 * C1;
        bigInt Z = B1 * D1;
        return multiple_by_ten_power(Y, 2*n) + multiple_by_ten_power(X - Y - Z, n) + Z;
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
