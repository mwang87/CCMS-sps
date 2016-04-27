#ifndef TWOVALUES_H
#define TWOVALUES_H

/**
* TODO: add description
*/
template<class T> class TwoValues {
public:
    /**
    * TODO: add description
    */
    T values[2];

    TwoValues() {
        values[0] = (T) 0;
        values[1] = (T) 0;
    }

    /**
    * TODO: add description
    */
    TwoValues(T a, T b) {
        values[0] = a;
        values[1] = b;
    }

    /**
    * TODO: add description
    *
    *@param other
    */
    TwoValues(const TwoValues<T> &other) {
        values[0] = other.values[0];
        values[1] = other.values[1];
    }

    /**
    * TODO: add description
    *
    *@param other
    *@return
    */
    TwoValues<T> &operator=(const TwoValues<T> &other) {
        values[0] = other.values[0];
        values[1] = other.values[1];
        return (*this);
    }

    /**
    * TODO: add description
    *
    *@param a
    *@param b
    *@return
    */
    TwoValues<T> &set(T a, T b) {
        values[0] = a;
        values[1] = b;
        return (*this);
    }

    /**
    * TODO: add description
    *
    *@param i
    *@return
    */
    T &operator[](int i) {
        return values[i];
    }

    const T &operator[](int i) const {
        return values[i];
    }

    /**
    * TODO: add description
    *
    *@param other
    *@return
    */
    bool operator<(const TwoValues<T> &other) {
        if (values[0] < other.values[0] or (values[0] == other.values[0]
                and values[1] < other.values[1]))
            return true;
        return false;
    }

    /**
    * TODO: add description
    *
    *@param other
    *@return
    */
    bool operator==(const TwoValues<T> &other) {
        return values[0] == other.values[0] and values[1] == other.values[1];
    }
    //    void operator<<(ostream &out) { out<<"("<<values[0]<<","<<values[1]<<")"; }
};

/**
* TODO: add description
*/
template<class T> bool operator<(const TwoValues<T> &a, const TwoValues<T> &b) {
    if (a.values[0] < b.values[0] or (a.values[0] == b.values[0]
            and a.values[1] < b.values[1]))
        return true;
    return false;
}

#endif
