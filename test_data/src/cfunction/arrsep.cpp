#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "mex.hpp"
#include "mexAdapter.hpp"


/*
    This function seperates the array for every "circular" distance between any two points are greater than gap. Make sure that the input range is 
    [-cirbnd, cirbnd]. If not, the distance between two points may be negative. (cirbnd = +Inf by default). For cirbnd < 0, then assume that there 
    is no circular bound, i.e. cirbnd = +Inf.

    Usage:
        out = arrsep(in, gap)
        out = arrsep(in, gap, cirbnd)
        out = arrsep(in, gap, cirbnd, tot)

        [input]:
            in:     the input array
            gap:    the minimum distance between adjacent point in the input array
            cirbnd: the circular bound for the array, this is used to measure the "distance" between angles, e.g. the "distance" between -0.49 and 0.49
                    is 0.1 + 0.1 = 0.2 with cirbnd = 0.5
            tot:    the small value to perform the comparsion between two float/bouble variables
        
        [output]:
            out:    the output array

        [default values]:
            cirbnd: inf
            tot:    1e-15

*/
// ref: https://www.mathworks.com/help/matlab/matlab_external/c-mex-source-file.html
// ref: https://www.mathworks.com/help/matlab/matlab_external/cpp-mex-api.html

using namespace matlab::data;
using matlab::mex::ArgumentList;


class MexFunction : public matlab::mex::Function
{
private:
    // Get array factory
    matlab::data::ArrayFactory factory;

    // Get pointer to engine
    const std::shared_ptr<matlab::engine::MATLABEngine>   matlabPtr = getEngine();

    // Create an output stream
    std::ostringstream stream;

    void displayOnMATLAB(std::ostringstream& stream)
    {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0, std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }

    /*
        DEBUG:
            stream << "in loop : " << i << std::endl;
            displayOnMATLAB(stream);
    */

   double gap;
   double cirbnd;
   double tot;

public:
    void operator()(ArgumentList outputs, ArgumentList inputs)
    {
        this->checkArguments(outputs, inputs);

        const size_t numRows = inputs[0].getDimensions()[0];
        const size_t numCols = inputs[0].getDimensions()[1];
        
        /* create data */
        const   matlab::data::TypedArray<double> arr_in     = std::move(inputs[0]);
                matlab::data::TypedArray<double> arr_out    = factory.createArray<double>({ numRows, numCols });

        gap = inputs[1][0];

        if (inputs.size() >= 3)
            cirbnd = inputs[2][0];
        else
            cirbnd = -1;

        if (inputs.size() >= 4)
            tot = inputs[3][0];
        else
            tot = 1e-15;

        // the following code seems to have some issue when input is double but get int
        // cirbnd = (inputs.size() >= 3) ? inputs[2][0] : -1;
        // thres  = (inputs.size() >= 4) ? inputs[3][0] : 1e-9;

        /* processing */
        this->separating(arr_in, arr_out);

        /* output */
        outputs[0] = std::move(arr_out);
    }

    void separating(const matlab::data::TypedArray<double>& arr_in, matlab::data::TypedArray<double>& arr_out)
    {     
        // ========== initialization ==========
        size_t N = arr_in.getNumberOfElements();

        // for simple cases:
        if (N == 1){
            arr_out[0] = arr_in[0];
            return;
        }
        else if ((cirbnd > 0) && (gap * N > 2 * cirbnd)){
            arr_out = std::move(arr_in);
            matlabPtr->feval(u"warning", 0, 
                std::vector<Array>({factory.createScalar("No feasible solution can be found.")}));
            return;
        }

        // ========== processing ==========
        std::vector<double>         new_arr(N);
        std::multimap<double, int>  sort_arr;
        bool    flag_too_close;
        int     close_idx;
        double  close_dst;
        
        new_arr[0] = arr_in[0];
        sort_arr.insert({arr_in[0], 0});

        for (int i = 1; i < N; i++){
            // insert to new array
            insertion_sort(new_arr, i+1, arr_in[i]);
            sort_arr.insert({arr_in[i], i});

            // check is there any points too close
            flag_too_close = true;

            while (flag_too_close){

                // 1. linear part
                if (find_too_close1(new_arr, i+1, close_idx, close_dst)){

                    // the "distance" needed to be inserted
                    close_dst = gap - close_dst;

                    // count contiguous left/right-hand side points
                    int cL = countL(new_arr, i+1, close_idx);
                    int cR = countR(new_arr, i+1, close_idx + 1);

                    // calculate the "push distance"
                    double pushL_dst = close_dst * cR / (cR + cL);
                    double pushR_dst = close_dst * cL / (cR + cL);

                    // deal with too small value
                    if(pushL_dst < tot) pushL_dst = tot;
                    if(pushR_dst < tot) pushR_dst = tot;
    
                    // push left/right
                    pushL(new_arr, i+1, close_idx, pushL_dst, cL);
                    pushR(new_arr, i+1, close_idx+1, pushR_dst, cR);
                }

                // 2. circular part
                else if (find_too_close2(new_arr, i+1, close_dst)){

                    // the "distance" needed to be inserted
                    close_dst = gap - close_dst;

                    // count contiguous left/right-hand side points
                    int cL = countL(new_arr, i+1, i);
                    int cR = countR(new_arr, i+1, 0);

                    // calculate the "push distance"
                    double pushL_dst = close_dst * cR / (cR + cL);
                    double pushR_dst = close_dst * cL / (cR + cL);

                    // deal with too small value
                    if(pushL_dst < tot) pushL_dst = tot;
                    if(pushR_dst < tot) pushR_dst = tot;

                    // push left/right
                    pushL(new_arr, i+1, i, pushL_dst, cL);      
                    pushR(new_arr, i+1, 0, pushR_dst, cR);

                }
                else{
                    flag_too_close = false;
                }
                
            }
        }

        /* if cirbnd > 0, deal with the circular boundary */
        if (cirbnd > 0){
            for (int i = 0; i < N; i++){
                if (new_arr[i] < -cirbnd){
                    new_arr[i] += 2 * cirbnd;
                }

                if (new_arr[i] > cirbnd){
                    new_arr[i] -= 2 * cirbnd;
                }
            }
        }

        // ========== move and unsort ==========
        int i = 0;
        for (auto& p : sort_arr){
            arr_out[p.second] = std::move(new_arr[i++]);
        }
    }

    // template<typename T>
    void circshift(std::vector<double>& v, const size_t N, int k) const
    {
        // Make sure k is in the range [0, N)
        k = k % N;
        if (k < 0) {
            // Convert negative k to a positive shift in the opposite direction
            k += N;
        }
        std::rotate(v.rbegin(), v.rbegin() + k, v.rend());
    }

    inline void pushL(std::vector<double>& arr, const size_t N, const int idx, const double pushL_dst, const int cL) const
    {
        for (int c = 0; c < cL; c++){
            arr[idx-c] -= pushL_dst;
        }
    }

    inline void pushR(std::vector<double>& arr, const size_t N, const int idx, const double pushR_dst, const int cR) const
    {
        for (int c = 0; c < cR; c++){
            arr[idx+c] += pushR_dst;
        }
    }

    inline int countL(const std::vector<double>& arr, const size_t N, int idx) const
    {
        /*
                    idx: 0 1 2 3 4 5
                    arr: 1 3 6 7 8 11
            count from:     <= ^      for idx = 3
            result: count_l = 2 (6 & 7)
        */
        int count_l = 1;

        while (idx > 0){
            if ((arr[idx] - arr[idx-1] - gap) < tot)
                count_l++;
            else
                return count_l;
            idx--;
        }
        return count_l;
    }

    inline int countR(const std::vector<double>& arr, const size_t N, int idx) const
    {
        /*
                    idx: 0 1 2 3 4 5
                    arr: 1 3 6 7 8 11
            count from:        ^ =>    for idx = 3
            result: count_r = 2 (7 & 8)
        */
        int count_r = 1;

        while (idx < N-1){
            if ((arr[idx+1] - arr[idx] - gap) < tot)
                count_r++;
            else
                return count_r;
            idx++;
        }
        return count_r;
    }

    bool find_too_close1(const std::vector<double>& arr, const size_t N, int &idx, double &val) const
    {
        // assume arr is in ascending order

        for (int k = 0; k < N-1; k++)
            if (arr[k+1] - arr[k] < gap){
                idx = k;
                val = arr[k+1] - arr[k];
                return true;
            }
        return false;
    }

    bool find_too_close2(const std::vector<double>& arr, const size_t N, double &val) const
    {
        // disable circular bound
        if (cirbnd <= 0)
            return false;

        // none or only one element
        if (N <= 1)
            return false;
        
        val = arr[N-1] - arr[0];
        val = std::min(val, 2 * cirbnd - val);

        return (val < gap);
    }

    void insertion_sort(std::vector<double>& arr, const int N, const double val)
    {
        // Insert the N-th element val to arr, it does not check if arr has the capacity of N
        // Note that N is one-index, i.e. N = 1, 2, 3, ...

        // find index of insertion
        int idx = 0;
        while (idx < N-1){
            if (val < arr[idx])
                break;
            idx++;
        }

        // insert value to arr
        for (int i = N-1; i > idx; i--){
            arr[i] = arr[i-1];
        }
        arr[idx] = val;
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs)
    {
        if (inputs.size() < 2){
            matlabPtr->feval(u"error", 0, 
                std::vector<Array>({factory.createScalar("At least two input required")}));
        }

        if (outputs.size() > 1){
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({factory.createScalar("Maximum one output supported")}));
        }

        if (inputs[0].getType() != ArrayType::DOUBLE || inputs[0].getType() == ArrayType::COMPLEX_DOUBLE ||
            inputs[1].getType() != ArrayType::DOUBLE || inputs[1].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({factory.createScalar("Input must be non-complex double elements")}));
        }

        // check first input dimension
        if (std::min(inputs[0].getDimensions()[0], inputs[0].getDimensions()[1]) != 1)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({factory.createScalar("Input must be an array")}));
        }

        // check second input dimension
        if (inputs[1].getNumberOfElements() != 1){
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({factory.createScalar("Second input must be a scalar")}));
        }

        // if there is third input
        if (inputs.size() > 2){
            // check type
            if (inputs[2].getType() != ArrayType::DOUBLE || inputs[2].getType() == ArrayType::COMPLEX_DOUBLE ||
                inputs[2].getType() != ArrayType::DOUBLE || inputs[2].getType() == ArrayType::COMPLEX_DOUBLE)
            {
                matlabPtr->feval(u"error", 0,
                    std::vector<Array>({factory.createScalar("Third input must be non-complex double elements")}));
            }
            // check dimension
            if (inputs[2].getNumberOfElements() != 1){
                matlabPtr->feval(u"error", 0,
                    std::vector<Array>({factory.createScalar("Third input must be a scalar")}));
            }
        }

        // if there is 4th input
        if (inputs.size() > 3){
            // check type
            if (inputs[3].getType() != ArrayType::DOUBLE || inputs[3].getType() == ArrayType::COMPLEX_DOUBLE ||
                inputs[3].getType() != ArrayType::DOUBLE || inputs[3].getType() == ArrayType::COMPLEX_DOUBLE)
            {
                matlabPtr->feval(u"error", 0,
                    std::vector<Array>({factory.createScalar("4th input must be non-complex double elements")}));
            }
            // check dimension
            if (inputs[3].getNumberOfElements() != 1){
                matlabPtr->feval(u"error", 0,
                    std::vector<Array>({factory.createScalar("4th input must be a scalar")}));
            }
        }
    }
};