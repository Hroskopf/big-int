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
            if(*this > other)
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

};
