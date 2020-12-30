#ifndef MDIMAGE_HPP_INCLUDED
#define MDIMAGE_HPP_INCLUDED

/** \brief Cimg.h: header-only C++ library for handling pictures
 *
 * Downloaded from https://framagit.org/dtschump/CImg/raw/master/CImg.h
 * it can save images as BMP picturs without requiring extra libraries
 */
#include "CImg/CImg.h"

/** \brief niftiio.hpp uses nifti1_io for NIfTI file I/O
 *         with a modern, C++ style interface
 */
#include "nifti/include/niftiio.hpp"

/** \brief mdimage_maths.hpp: maths often used in images
 *
 * Including Gaussian filters/noise and 3D/4D vectors
 */
#include "mdimage_maths.hpp"


#include <queue>
#include <random>
#include <limits>
#include <numeric>
#include <algorithm>
#include <initializer_list>

/** \brief NO/DO read data enumeration for easily readable code
 */
enum readdata {
    NO_READ_DATA,
    DO_READ_DATA
};



/** \brief mdi: namespace for multidimensional images
 *
 *  This namespace contains the mdimage class for
 *  creating and processing multidimensional images.
 */
namespace mdi {




/** \brief filterline for 1-dimensional convolution
 *
 * signal and filter are both of type std::vector
 */
template <typename T, typename U>
void filterline ( std::vector<T>& signal, const std::vector<U>& filter ) {

    // do circular for now -> after signal[size -1] comes signal[0]
    vector <T> workspace ( signal.size() );
    size_t flen2 = filter.size() / 2;

    for ( size_t si = 0; si < signal.size(); si++ )
        for ( size_t fi = 0; fi < filter.size(); fi++ )
            workspace [ (si + flen2) % signal.size() ] += signal [ (si + fi) % signal.size() ] * filter [ fi ];

    std::copy ( workspace.begin(), workspace.end(), signal.begin() );

}

/** \brief fwtline for 1-dimensional fast wavelet transform
 *
 * signal and filterh are both of type std::vector
 * level is of type size_t
 */
template <typename T, typename U>
void fwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level ) {

    // do circular for now -> after signal[size -1] comes signal[0]
    size_t
    newlength = signal.size();
    std::vector <U>
    filterg = ( filterh );

    // filter g not only the reverse of h, every second coefficient must be - as well
    std::reverse ( filterg.begin(), filterg.end() );
    for ( size_t i=1; i < filterg.size(); i+=2 )
        filterg[i] *= -1;

    for ( size_t l=0; l<level; l++ ) {

        std::vector <T>
        tmp ( newlength, 0 );

        newlength /= 2;

        for ( size_t wi=0; wi<newlength; wi++ )
            for ( size_t fi=0; fi<filterh.size(); fi++ ) {

                tmp[wi            ] += signal [ (2*wi + fi) % signal.size() ] * filterh [fi];
                tmp[wi + newlength] += signal [ (2*wi + fi) % signal.size() ] * filterg [fi];

            } // for wi

        std::copy ( tmp.begin(), tmp.end(), signal.begin() );

    } // for l

} // fwtline

/** \brief fwtline for 1-dimensional fast inverse wavelet transform (reconstruction)
 *
 * signal and filterh are both of type std::vector
 * level is of type size_t
 */
template <typename T, typename U>
void ifwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level ) {

    // do circular for now -> after signal[size -1] comes signal[0]
    size_t
    newlength = signal.size();
    std::vector <U>
    filterg = ( filterh );

    // filter g not only the reverse of h, every second coefficient must be 0 as well
    std::reverse ( filterg.begin(), filterg.end() );
    for ( size_t i=1; i < filterg.size(); i+=2 )
        filterg[i] *= -1;

    // length of approximation is 2^level times smaller
    for ( size_t i=0; i<level; i++ )
        newlength /= 2;

    for ( size_t l=0; l<level; l++ ) {

        std::vector <T>
        tmpc ( newlength*2, 0 );
        std::vector <T>
        tmpd ( newlength*2, 0 );

        // upsample c and d from current level
        for ( size_t i=0; i<newlength; i++ ) {
            tmpc [ 2 * i ] = signal [ i             ];
            tmpd [ 2 * i ] = signal [ newlength + i ];
        }

        std::fill( signal.begin(), signal.end(), 0 );
        newlength *=2;

        for ( size_t wi=0; wi<newlength; wi++ )
            for ( size_t fi=0; fi<filterh.size(); fi++ ) {

                signal[wi] += tmpc [ (wi - fi) % tmpc.size() ] * filterh[fi];
                signal[wi] += tmpd [ (wi - fi) % tmpd.size() ] * filterg[fi];

            } // for wi

    } // for l

} // ifwtline















/** \brief class template imageNd for n-dimensional images
 *
 * voxel type T can be any scalar type (char, short, int, float, double)
 */
template <typename T>
class imageNd {

    
    
  public:

    /** \brief use value_type for pixels
     *
     * To be able to use value_type (this name is also
     * used in STL containers) we make it public first
     */
    typedef
    T value_type;



  private:

    /** \brief use value_type for pixels
     *
     * data:    vector with the pixel data
     * sizes:   vector with the dimensions
     * strides: vector with the strides through
     *          the data in each dimension
     * header:  pointer to a NIfTI header
     *          for reading / writing files
     *
     */
    std::vector <size_t>
    sizes;

    std::vector <size_t>
    strides;

    nifti_image
    *header = NULL;

    std::vector <value_type>
    data;



  public:

    /** \brief default constructor
     */
    imageNd() {}

    /** \brief destructor
     *
     * destructor: clears vectors and (if
     * required) frees pointer to header
     */
    ~imageNd() {
        
        data.resize    (0);
        sizes.resize   (0);
        strides.resize (0);
        if ( header != NULL )
            free (header);

    }

    /** \brief constructor for an empty image
     *                            with sizes
     *
     * constructor for an empty image (no
     * nifti header) with given dimensions
     */
    imageNd ( std::initializer_list <size_t> dims ) {

        sizes.resize   ( dims.size()   );
        strides.resize ( dims.size()+1 );

        for ( size_t i=1; i<=sizes.size(); i++ ) {
            sizes [i-1] = dims.begin()[i];
            strides [i] = strides[i-1] * sizes[i-1];
        }

        data.resize ( *strides.rbegin() );

    }

    /** \brief (deep) copy constructor
     *
     * copies data, sizes and header from en existing imageNd
     */
    imageNd ( const imageNd& rhs ) {

        if (rhs.header != NULL)
            header = nifti_copy_nim_info (rhs.header);

        data.resize( rhs.data.size() );
        sizes.resize( rhs.sizes.size() );
        strides.resize( rhs.strides.size() );

        std::copy ( rhs.data.begin(),    rhs.data.end(),    data.begin()    );
        std::copy ( rhs.sizes.begin(),   rhs.sizes.end(),   sizes.begin()   );
        std::copy ( rhs.strides.begin(), rhs.strides.end(), strides.begin() );

    }

    /** \brief assignment operator
     *
     * assigns data, sizes and header of the
     * right-hand side (RHS) image to (*this)
     */
    const imageNd<T>& operator= ( imageNd <T>& rhs ) {

        if ( this != &rhs ) {

            // just to make sure we don't leave stuff
            if ( this->header != NULL ) {
                free ( header );
                header = NULL;
            }

            if (rhs.header != NULL)
                header = nifti_copy_nim_info( rhs.header );

            // just to make sure we don't leave stuff
            data.resize    ( rhs.data.size()    );
            sizes.resize   ( rhs.sizes.size()   );
            strides.resize ( rhs.strides.size() );

            std::copy ( rhs.data.begin(),    rhs.data.end(),    data.begin()    );
            std::copy ( rhs.sizes.begin(),   rhs.sizes.end(),   sizes.begin()   );
            std::copy ( rhs.strides.begin(), rhs.strides.end(), strides.begin() );

        } // if this != rhs

        return *this;

    } // assignment

    /** \brief constructor from a NIfTI image
     *
     * This constructor reads an n-dimensional image
     * ( for NIfTI 1 <= n <= 7 ) from a NIfTI file and
     * assigns its data, sizes and header to (*this).
     */
    imageNd ( std::string filename ) {

        header = nifti_image_read ( filename.c_str(),
                                    DO_READ_DATA );

        sizes.resize   ( header -> dim[0]   );
        strides.resize ( header -> dim[0]+1 );

        // make the array 'strides' so that it uses the last dimension as well
        strides [0] = 1;
        for ( size_t i=1; i<=sizes.size(); i++ ) {
            sizes [i-1] = header -> dim [i];
            strides [i] = strides[i-1] * sizes[i-1];
        }

        data.resize ( *(strides.rbegin()) ); // the end of strides holds the image's size

        aip::getNiftiBricks ( header,
                              header -> data,
                              data.size(),
                              &data );

    }

    /** \brief change the data type for writing a NIfTI file
     *
     * The data type is changed to e.g. float for
     * more flexibility (at the cost of file size).
     */
    void setNIIdatatype ( unsigned dtype ) {
        
        header-> datatype = dtype;
        nifti_datatype_sizes( header->datatype, &header->nbyper, &header->swapsize ) ;
        
    }

    /** \brief write to a NIfTI image
     *
     * Writes the contents of (*this) to a NIfTI
     * file using the data and header information.
     * The data type is changed to float for more
     * flexibility (at the cost of file size).
     *
     * NIfTI routines are run on a temporary copy
     *       as they seem to be memory-unsafe
     */
    void saveNII ( std::string filename ) {
        
        nifti_set_filenames ( header, filename.c_str(), 0, 0 );
        setNiftiBricks      ( header,            &data       );
        nifti_image_write( header );
        
    }



    /** \brief operator() for positional addressing
     *
     * for an existing image In, this operator can
     * be used with multidimensional co-ordinates to
     * indicate position, so In({x,y,z}) instead of
     * In.data[x + y*sizes[1] + z*sizes[1]*sizes[2]].
     *
     * This operator is for reading only.
     */
    value_type const operator() ( std::initializer_list < size_t > const& indices ) const {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ];
    }

    /** \brief operator() for positional addressing
     *
     * This operator is for modifying the data.
     */
    value_type& operator() ( std::initializer_list < size_t > const& indices ) {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ] ;
    }

    /** \brief operator[] for positional addressing
     *
     * for an existing image In, this operator can
     * be used with multidimensional co-ordinates to
     * indicate position, so In[{x,y,z}] instead of
     * In.data[x + y*sizes[1] + z*sizes[1]*sizes[2]].
     *
     * This operator is for reading only.
     */
    value_type const operator[] ( std::initializer_list < size_t > const& indices ) const {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ];
    }

    /** \brief operator[] for positional addressing
     *
     * This operator is for modifying the data.
     */
    value_type& operator[] ( std::initializer_list < size_t > const& indices ) {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ] ;
    }

    /** brief compute indices at offset
     *
     * inverse of positional addressing - given a position
     * in the 1D data vector, what are its multidimensional
     * indices?
     */
    const std::vector <size_t> indices_at_offset ( size_t pos1d ) {

        auto p = pos1d;
        std::vector <size_t> out ( sizes.size() );

        for ( size_t d = sizes.size()-1; d>0; d-- ) {
            out[d]  =      p / strides[d];
            p      -= out[d] * strides[d];
        }
        out[0] = p;

        return out;

    }


    /** \brief operators += for scalar and imageNd, repectively
     *
     */
    const imageNd<T>& operator+= ( const value_type& rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] += rhs;
        return (*this);
    }
    template <typename U>
    const imageNd<T>& operator+= ( const imageNd<U>& rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] += rhs.data[s];
        return (*this);
    }

    /** \brief operator + for templated types
     *
     * in this case, types for which += has been defined
     */
    template <typename U>
    const imageNd<T> operator+ ( const U& rhs ) {
        imageNd out(*this);
        out += rhs;
        return out;
    }



    /** \brief operators *= for scalar and imageNd, repectively
     *
     */
    const imageNd<T>& operator*= ( const value_type& rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] *= rhs;
        return (*this);
    }
    template <typename U>
    const imageNd<T>& operator*= ( const imageNd<U>& rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] *= rhs.data[s];
        return (*this);
    }

    /** \brief operator * for templated types
     *
     * in this case, types for which *= has been defined
     */
    template <typename U>
    const imageNd<T> operator* ( const U& rhs ) {
        imageNd out(*this);
        out *= rhs;
        return out;
    }



    /** \brief operators -= for scalar and imageNd, repectively
     *
     */
    const imageNd<T>& operator-= ( const value_type rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] -= rhs;
        return (*this);
    }
    template <typename U>
    const imageNd<T>& operator-= ( const imageNd<U>& rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] -= rhs.data[s];
        return (*this);
    }

    /** \brief operator - for templated types
     *
     * in this case, types for which -= has been defined
     */
    template <typename U>
    const imageNd<T> operator- ( const U& rhs ) {
        imageNd out(*this);
        out -= rhs;
        return out;
    }

    /* \brief operator - without operand
     *
     * negate yourself
     */
    const imageNd<T> operator- ( void ) {
        imageNd out(*this);
        out *= -1;
        return out;
    }


    /** \brief operators /= for scalar and imageNd, repectively
     *
     */
    const imageNd<T>& operator/= ( const value_type rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] /= rhs;
        return (*this);
    }
    template <typename U>
    const imageNd<T>& operator/= ( const imageNd<U>& rhs ) {
        for (size_t s = 0; s < data.size(); s++ )
            data[s] /= rhs.data[s];
        return (*this);
    }

    /** \brief operator / for templated types
     *
     * in this case, types for which /= has been defined
     */
    template <typename U>
    const imageNd<T> operator/ ( const U& rhs ) {
        imageNd<T> out(*this);
        out /= rhs;
        return out;
    }



    /** \brief reciprocal function, returns 1/c
     *
     * Return an image with 1/c for all coefficients c.
     * This function is used for dividing by an image.
     */
    const imageNd<T>& reciprocal ( const imageNd<T>& rhs ) {
        imageNd<T> out(*rhs);
        for (size_t s = 0; s < out.data.size(); s++ )
            out.data[s] = 1 / out.data[s];
        return out;
    }



    /** \brief getsize() - returns the dimensions
     *
     * return dimensions as a std::vector <size_t>
     */
    std::vector<size_t> getsize (          ) {
        return sizes;
    }



    /** \brief getsize() - returns one dimension
     *
     * return dimension given its index in sizes
     */
    size_t              getsize ( size_t s ) {
        return sizes[s];
    }



    /** \brief getdatasize() - returns the number of intensities
     *
     * the size of the vector 'data'
     */
    size_t              getdatasize ( ) {
        return std::accumulate( sizes.begin(), sizes.end(), 1, std::multiplies<size_t>() );
    }



    /** \brief getdata_ptr() - returns the address of the data vector
     *
     * this is a pointer to a vector -- use with care
     */
    std::vector <value_type>*          getdata_ptr ( ) {
        return &data;
    }



    /** \brief getdata_array() - returns the address of the 1st data value
     *
     * this is a pointer to an array -- use with care
     */
    value_type*                        getdata_array ( ) {
        return data.data();
    }

    void setdata_array ( value_type* newdata ) {
        for ( size_t i = 0; i < data.size(); i++ )
            data[i] = newdata[i];
    }





    /** \brief reshape() - change dimensions
     *
     * Changes the sizes and strides vectors.
     * Only works if total #elements does not change.
     */
    void reshape ( std::initializer_list < size_t > const& newsizes ) {

        if ( std::accumulate( newsizes.begin(), newsizes.end(), 1, std::multiplies<size_t>() ) ==
                std::accumulate(    sizes.begin(),    sizes.end(), 1, std::multiplies<size_t>() )
           ) {

            sizes=newsizes;
            for ( size_t i=1; i<=sizes.size(); i++ ) {
                strides [i] =  strides[i-1] * sizes[i-1];
            }

        } else

            std::cout << "reshape impossible because it would change image size";

    } // reshape



    /** \brief addNormalNoise() add normally distributed noise
     *
     * Uses standard library functions to produce random numbers
     * Parameters mu and sigma are doubles, used in the expected way
     */
    void addNormalNoise ( double mu, double sigma ) {

        // random device class instance, source of 'true' randomness for initializing random seed
        std::random_device randomdevice{};

        // Mersenne twister PRNG, initialized with seed from previous random device instance
        std::mt19937 engine { randomdevice() };

        // set up the distribution
        std::normal_distribution <double> normaldist ( mu, sigma );

        // add noise to the data
        // N (mu, sigma) all with <double> values
        // (non-float leads to unexpected behaviour)
        for ( size_t i = 0; i < data.size(); i++ )
            data[i] += normaldist ( engine );

    }

    /** \brief addRicianNoise() add Rician distributed noise
     *
     * Uses standard library functions to produce random numbers
     * Parameters mu and sigma are doubles, used in the expected way
     * for two normal distributions.
     */
    void addRicianNoise ( double mu, double sigma ) {

        // random device class instance, source of 'true' randomness for initializing random seed
        std::random_device randomdevice{};

        // Mersenne twister PRNG, initialized with seed from previous random device instance
        std::mt19937 engine { randomdevice() };

        // set up the distribution
        std::normal_distribution <double> normaldist ( mu, sigma );

        // add Rician noise to the data using 2 normally distributed noise values
        for ( size_t i = 0; i < data.size(); i++ ) {

            double n1 = data[i] + normaldist ( engine );
            double n2           = normaldist ( engine );
            data [i] = sqrt ( ( n1 * n1 ) + ( n2 * n2 ) );

        }

    }



    /** \brief get3Dline() returns a 3D line (along a dimension) from a volume
     *
     * Result is a std::vector of value_type
     *
     * The line is sampled along dimension dim ( 0, 1 or 2 respectively)
     * at position pos1, pos2 in the other dimensions
     * ( 1 and 2, 0 and 2, or 0 and 1, respectively).
     */
    const std::vector <value_type> get3Dline ( size_t dim, size_t pos1, size_t pos2,
            size_t linestart = 0, size_t lineend = UINT_MAX ) {

        std::vector <value_type> out;

        if (sizes.size() != 3) {

            std::cout << "3D lines must be selected from a 3D image\n";

        } else {

            size_t
            step    = strides[dim],
            slicex  = (dim>0) ? 0 : 1,
            slicey  = (dim>1) ? 1 : 2,
            line0   = std::max <size_t> ( linestart, 0          ),
            line1   = std::min <size_t> ( lineend,   sizes[dim] );

            value_type*
            dptr    = data.data() + pos1 * strides[slicex] + pos2 * strides[slicey];

            out.resize ( line1-line0 );

            for ( size_t i = line0; i < line1; i++, dptr+=step )
                out[i] = *dptr;

        } // if sizes

        return out;

    } // get3Dline



    /** \brief get3Dline() puts a 3D line (along a dimension) in a volume
     *
     * Input is a std::vector of value_type
     *
     * The line is inserted along dimension dim ( 0, 1 or 2 respectively)
     * at position pos1, pos2 in the other dimensions
     * ( 1 and 2, 0 and 2, or 0 and 1, respectively).
     */
    void set3Dline ( std::vector <value_type>& in,
                     size_t dim, size_t pos1, size_t pos2,
                     size_t linestart = 0, size_t lineend = UINT_MAX ) {

        if (sizes.size() != 3) {

            std::cout << "3D lines must be selected from a 3D image\n";

        } else {

            size_t
            step = strides[dim],
            slicex = (dim>0) ? 0 : 1,
            slicey = (dim>1) ? 1 : 2,
            line0 = std::max <size_t> ( linestart, 0          ),
            line1 = std::min <size_t> ( line0 + lineend, sizes[dim] );

            value_type* dptr = data.data() + pos1 * strides[slicex] + pos2 * strides[slicey];

            for ( size_t i = line0; i < line1; i++, dptr+=step )
                *dptr = in[i];

        } // if sizes

    } // set3Dline



    /** \brief getSlice() returns a 2D slice from a volume
     *
     * Result is an imageNd of value_type
     * which is a copy of slice no. <sli>
     * along dimension <dim> (0, 1, or 2)
     * and at position <sli>.
     */
    const imageNd<value_type> getSlice( size_t dim, size_t sli, std::string filename="" ) {
        // get a slice from a volume
        // optional: write it out as a .bmp file
        //

        imageNd<value_type> out;

        if (sizes.size() != 3) {

            std::cout << "slices can only be selected from a 3D image\n";

        } else {

            // slice sizes are called slicex (lowest 3D dim index) and slicey (highest)
            size_t slicex = (dim>0) ? 0 : 1;
            size_t slicey = (dim>1) ? 1 : 2;

            // set sizes for sizes, strides and data start with 3D
            out.sizes     = {    sizes[slicex], sizes[slicey], 1                                             };
            out.strides   = { 1, sizes[slicex], sizes[slicex] * sizes[slicey], sizes[slicex] * sizes[slicey] };
            out.data.resize (    *out.strides.rbegin()                                                       );

            // fill the slice by calling get3Dline (from folume) for each y line in the slice
            // loop over highest (outer == slower with largest strides) dimension first
            //
            // dim x -> yz slice, slicex 1, slicey 2 -> lines in y, loop over z -> line pos [ sli z ] = [  sli ypos ]
            // dim y -> xz slice, slicex 0, slicey 2 -> lines in x, loop over z -> line pos [ sli z ] = [  sli ypos ]
            // dim z -> xy slice, slicex 0, slicey 1 -> lines in x, loop over y -> line pos [ y sli ] = [ ypos  sli ]
            for ( size_t ypos=0; ypos<sizes[slicey]; ypos++ ) {

                // position where the line is taken:
                //
                // an x line is taken from an y,z position
                size_t linx = ( slicey>1 ) ?  sli : ypos;
                size_t liny = ( slicey>1 ) ? ypos :  sli;

                std::vector<value_type> sli_line = get3Dline( slicex, linx, liny );
                out.set3Dline ( sli_line, 0, ypos, 0); // x line (0), put in y position liny, 'z' position 0

            } // for ypos

        } // if sizes

        if ( !filename.empty() ) {

            cimg_library::CImg<value_type>*
            my_bitmap = new cimg_library::CImg<value_type> (out.data.data(),
                    out.sizes[0],
                    out.sizes[1],
                    1, 1, true);
            my_bitmap->rotate(180);
            my_bitmap->save_bmp( filename.c_str() );
            delete ( my_bitmap );

        } // if filename

        out.reshape( { out.sizes[0], out.sizes[1] } ); // remove dimension 3 (which is 1) of the output slice
        return out;

    } // getSlice

    /** \brief setSlice() insert a 2D slice into a volume
     *
     * Slice no. <sli> along
     * dimension <dim> (0, 1, or 2)
     * of the current object is copied from the input
     */
    void setSlice ( imageNd<value_type>& input, size_t dim, size_t sli ) {

        if (sizes.size() != 3) {

            std::cout << "slices can only be inserted into a 3D image\n";

        } else {

            // slice sizes are called slicex (lowest 3D dim index) and slicey (highest)
            size_t slicex = (dim>0) ? 0 : 1;
            size_t slicey = (dim>1) ? 1 : 2;

            // check if input sizes match the current imagNd dimensions
            if ( ( input.sizes[0] != sizes[slicex] ) | ( input.sizes[1] != sizes[slicey] ) ) {

                std::cout << "input slice deimensions do not match volume slice size \n";

            }

            // briefly make our slice 3D for using get3Dline()
            input.reshape ( input.sizes[0], input.size[1], 1 );

            for ( size_t ypos=0; ypos<sizes[slicey]; ypos++ ) {

                size_t linx = ( slicey>1 ) ?  sli : ypos;
                size_t liny = ( slicey>1 ) ? ypos :  sli;

                std::vector<value_type> sli_line = input.get3Dline ( sli_line, 0, ypos, 0);
                set3Dline( slicex, linx, liny );

            } // for ypos

            // make our slice 2D again
            input.reshape ( input.sizes[0], input.size[1] );

        } // if sizes

    } // getSlice



    /** \brief getsubvolume() returns a 3D subvolume from a 3D image
     *
     * Ranges are given as xmin, xmax, ymin, ymax, zmin, zmax
     * Input and output are both of type imageNd, and this
     * routine works exclusively with 3D data
     */
    const imageNd <value_type> getsubvolume (size_t startx, size_t endx,
            size_t starty, size_t endy,
            size_t startz, size_t endz) {

        imageNd <value_type> out;

        if (sizes.size() != 3) {

            std::cout << "subvolumes can only be selected from a 3D image\n";

        } else {

            size_t  x0 = std::max<size_t> ( startx, 0 ),
                    y0 = std::max<size_t> ( starty, 0 ),
                    z0 = std::max<size_t> ( startz, 0 ),
                    x1 = std::min<size_t> ( endx, sizes[0] ),
                    y1 = std::min<size_t> ( endy, sizes[1] ),
                    z1 = std::min<size_t> ( endz, sizes[2] );

            out.sizes   = { std::max<size_t> (x1 - x0, 1),
                            std::max<size_t> (y1 - y0, 1),
                            std::max<size_t> (z1 - z0, 1)
                          };
            out.strides = { 1, out.sizes[0], out.sizes[1] * out.sizes[0], out.sizes[2] * out.sizes [1] * out.sizes[0] };
            out.data.resize( *out.strides.rbegin() );

            value_type *dptr = out.data.data();

            for ( size_t z=z0; z<z1; z++ )
                for ( size_t y=y0; y<y1; y++ )
                    for ( size_t x=x0; x<x1; x++ )
                        *dptr++ = operator[] ( { x, y, z } );

        } // if sizes

        return out;

    } // getsubvolume

    /** \brief setsubvolume() inserts a 3D subvolume into a 3D image
     *
     * Ranges are given as xmin, ymin, zmin for where to insert
     * Source and destination are both of type imageNd, and this
     * routine works exclusively with 3D data
     */
    void setsubvolume ( imageNd <value_type>& in,
                        size_t startx,
                        size_t starty,
                        size_t startz ) {

        if ( (sizes.size() != 3) | (in.sizes.size() != 3) ) {

            std::cout << "only 3D can be put in only 3D images\n";

        } else {

            size_t  x0 = std::max<size_t> ( startx, 0 ),
                    y0 = std::max<size_t> ( starty, 0 ),
                    z0 = std::max<size_t> ( startz, 0 ),
                    x1 = std::min<size_t> ( startx + in.sizes[0], sizes[0] ),
                    y1 = std::min<size_t> ( starty + in.sizes[1], sizes[1] ),
                    z1 = std::min<size_t> ( startz + in.sizes[2], sizes[2] );

            value_type *dptr = &in.data[0];

            for ( size_t z=z0; z<z1; z++ )
                for ( size_t y=y0; y<y1; y++ )
                    for ( size_t x=x0; x<x1; x++ )
                        operator[] ({ x, y, z }) = *dptr++;

        } // if sizes

    } // getsubvolume



    /** \brief filter() filter along dimension {0, 1 or 2} in a 3D image
     *
     * The filter is given as a numrical std::vector
     * this method uses the function filterline
     * (outside this class)
     */
    void filter ( std::vector<double> filt, size_t dim ) {

        if (sizes.size() != 3) {

            std::cout << "currently filter works only for 3D images\n";

        } else {

            size_t slicex = (dim>0) ? 0 : 1;
            size_t slicey = (dim>1) ? 1 : 2;

            for ( size_t posx = 0; posx < sizes[slicex]; posx++ )

                for ( size_t posy = 0; posy < sizes[slicey]; posy++ ) {

                    std::vector <value_type> sign = get3Dline( dim, posx, posy );
                    filterline ( sign, filt );
                    set3Dline( sign, dim, posx, posy);

                } // for posy


        } // if sizes

    } //filter

    /** \brief fwt_slice() apply the fwt to the 2 dimensions of each slice
     *
     * Slice direction is taken as the 1st dimension (this is not standard!).
     * Dimensions 0 and 1 in each slice are dimension 1 and 2 of the input.
     * This methos uses fwtline (outside this class), wavelet basis is given
     * by its name ("Daubecchies4", "Haar").
     */
    void fwt_slices ( std::string wname, size_t level ) {

        // get filter as a vector not a name
        std::vector <double> wvec = owfilter( wname );

        // set sizes
        size_t xmax = sizes[0];
        size_t ymax = sizes[1];
        size_t zmax = sizes[2];

        // loop over levels
        for ( size_t l=0; l<level; l++ ) {

            // loop over x slices, do fwt in z then y
            for ( size_t x=0; x<xmax; x++ ) {

                for ( size_t y=0; y<ymax; y++ ) {
                    std::vector <value_type> dwt1 = get3Dline( 2, x, y, 0, zmax );
                    fwtline( dwt1, wvec, 1 );
                    set3Dline( dwt1, 2, x, y, 0, zmax );
                } // for y

                for ( size_t z=0; z<zmax; z++ ) {
                    std::vector <value_type> dwt1 = get3Dline( 1, x, z, 0, ymax );
                    fwtline( dwt1, wvec, 1 );
                    set3Dline( dwt1, 1, x, z, 0, ymax );
                } // for y

            } // for x

            ymax /= 2;
            zmax /= 2;

        } // for l

    } // fwtslices

    /** \brief ifwt_slice() apply the fwt to the 2 dimensions of each slice
     *
     * Slice direction is taken as the 1st dimension (this is not standard!).
     * Dimensions 0 and 1 in each slice are dimension 1 and 2 of the input.
     * This methos uses ifwtline (outside this class), wavelet basis is given
     * by its name ("Daubecchies4", "Haar").
     */
    void ifwt_slices ( std::string wname, size_t level ) {

        // get filter as a vector not a name
        std::vector <double> wvec = owfilter( wname );

        size_t fac=1;
        for (size_t i=1; i<level; i++)
            fac *= 2;

        // set sizes
        size_t xmax = sizes[0];
        size_t ymax = sizes[1]/fac;
        size_t zmax = sizes[2]/fac;

        // loop over levels
        for ( size_t l=0; l<level; l++ ) {

            // loop over x slices, do fwt in z then y
            for ( size_t x=0; x<xmax; x++ ) {

                for ( size_t y=0; y<ymax; y++ ) {
                    std::vector <value_type> dwt1 = get3Dline( 2, x, y, 0, zmax );
                    ifwtline( dwt1, wvec, 1 );
                    set3Dline( dwt1, 2, x, y, 0, zmax );
                } // for y

                for ( size_t z=0; z<zmax; z++ ) {
                    std::vector <value_type> dwt1 = get3Dline( 1, x, z, 0, ymax );
                    ifwtline( dwt1, wvec, 1 );
                    set3Dline( dwt1, 1, x, z, 0, ymax );
                } // for y

            } // for x

            ymax *= 2;
            zmax *= 2;

        } // for l

    } // ifwtslices

    /** \brief waveletVisuThresh() apply hard VisuThresh shrinkage
     *         to 2D wavelet channels of each sllice in a 3D image
     *
     * This method is only useful on images after applying
     * the fwt_slices() method, see above, at the same level
     */
    template<typename V>
    void waveletVisuThresh ( size_t level, V booster ) {

        // for 2D wavelet channels in slices, apply the VisuThresh
        // denoising threshold to the detail coefficients. This is
        // defined as sqrt(2 * log(|c|)) wher c is the size of the
        // wavelet channel. These channels are the same size in all
        // slices, so we can select them together.

        size_t xmax = sizes[0];
        size_t ymax = sizes[1];
        size_t zmax = sizes[2];

        for ( size_t l = 0; l< level; l++ ) {

            ymax /= 2;
            zmax /= 2;

            double threshold = booster * sqrt ( 2 * log (zmax) / log (2) );

            {
                // y only
                imageNd <value_type> workspace = getsubvolume(  0, xmax,
                                                 ymax, 2*ymax,
                                                 0, zmax );
                for ( size_t i=0; i<workspace.data.size(); i++)
                    workspace.data[i] *= (workspace.data[i]>threshold);
                setsubvolume( workspace, 0, ymax, 0 );
            }

            {
                // z only
                imageNd <value_type> workspace = getsubvolume(  0, xmax,
                                                 0, ymax,
                                                 zmax, 2*zmax );
                for ( size_t i=0; i<workspace.data.size(); i++)
                    workspace.data[i] *= (workspace.data[i]>threshold);
                setsubvolume( workspace, 0, 0, zmax );
            }

            {
                // both y and z
                imageNd <value_type> workspace = getsubvolume(  0, xmax,
                                                 ymax, 2*ymax,
                                                 zmax, 2*zmax );
                for ( size_t i=0; i<workspace.data.size(); i++)
                    workspace.data[i] *= (workspace.data[i]>threshold);
                setsubvolume( workspace, 0, ymax, zmax );
            }

        } // for l

    } // waveletVisuThresh















    int GetNNeigbours(int ip, int* NeighbourIndices, int ndim, size_t* dimensions) {
        if(ndim<=0 || ndim>3) {
            std::cout<<"ERROR: MImage::GetNNeigbours(). ndim (="<< ndim <<") out of range. \n";
            return 0;
        }
        if(NeighbourIndices==NULL) {
            std::cout<<"ERROR: MImage::GetNNeigbours(). Invalid NULL argument. \n";
            return 0;
        }

        int dimx  = dimensions[0];
        int dimy  = dimensions[1];
        int dimz  = dimensions[2];
        int dimxy = dimx*dimy;

// Test for out of range
        int ix = ndim>0 ? (  ip    % dimx       ) : 0;
        int iy = ndim>1 ? (((ip-ix)/ dimx)%dimy ) : 0;
        int iz = ndim>2 ? (  ip    / dimxy      ) : 0;

        if((ndim>0 && (ix<0 || ix>=dimx))  ||
                (ndim>1 && (iy<0 || iy>=dimy))  ||
                (ndim>2 && (iz<0 || iz>=dimz))) {
            std::cout<<"ERROR: MImage::GetNNeigbours(). point index out of range (ix, iy, iz) = ("<< ix <<", "<< iy <<" "<< iz << ")\n";
            return 0;
        }

        int NNeig = 0;
        if(ndim>0 && dimx>1) {
            if(ix>0     )
                NeighbourIndices[NNeig++]=ip-1;
            if(ix<dimx-1)
                NeighbourIndices[NNeig++]=ip+1;
        }
        if(ndim>1 && dimy>1) {
            if(iy>0     )
                NeighbourIndices[NNeig++]=ip-dimx;
            if(iy<dimy-1)
                NeighbourIndices[NNeig++]=ip+dimx;
        }
        if(ndim>2 && dimz>1) {
            if(iz>0     )
                NeighbourIndices[NNeig++]=ip-dimxy;
            if(iz<dimz-1)
                NeighbourIndices[NNeig++]=ip+dimxy;
        }
        return NNeig;
    }


















    bool GetWatershedImage() {

        if ( sizes.size() > 3 ) {
            std::cout<<"ERROR: MImage::GetWatershedImage(). Invalid dimensionality ("<< sizes.size() <<"). \n";
            return false;
        }

        value_type max_value=*std::max_element(data.begin(),data.end());
        value_type min_value=*std::min_element(data.begin(),data.end());
        size_t ndim=sizes.size();
        size_t NP=getdatasize();

        std::vector<int>
        Index   (NP),
                Dist    (NP),
                Label   (NP),
                Hist    (max_value + 2),
                CHist   (max_value + 2),
                NeigArr (200);

        // check if image needs inverting and do so if yes
        // (if pixel [0] has lower value than maximum/2 ?)
        if ( data[0] < (max_value/2) ) {
            std::cout << "inverting ... \n";
            for ( size_t i=0; i< NP; i++ )
                data [i] = max_value - min_value - data[i];
        }

        // build the histogram
        for (unsigned n=0; n < NP; n++)
            Hist[data[n]]++;

        // build the cumulative histogram (differs from histogram after index 0)
        for (unsigned k=1; k < max_value+1; k++)
            CHist[k] = CHist[k-1] + Hist[k-1];

        // label point based on value in cumulative histogram -- increasing index to number within intensity
        for (unsigned n=0; n < NP; n++)
            Index[CHist[data[n]]++] = n;

        // subtract histogram from cumulative after labelling
        for (unsigned k=0; k< max_value+1; k++)
            CHist[k] -= Hist[k]; // restore cumulative histogram

        CHist[max_value+1] = NP; // this was still 0

        const int LABELINIT  =   -1;
        const int MASK       =   -2;
        const int WSHED      =    0;
        const int FICTITIOUS =   -3;

        // initialise labels
        for ( unsigned n=0; n< NP; n++)
            Label[n] = LABELINIT;

        std::queue<int> fifoQueue;
        int curlab = 0;

        // Geodesic SKIZ of level h-1 inside level h. INCLUDE LAST LEVEL!
        for( value_type h = min_value; h<=max_value; h++) {
            for( int pixelIndex = CHist[h]; pixelIndex < CHist[h+1]; pixelIndex++) { //mask all pixels at level h
                int   ip  = Index[pixelIndex];
                Label[ip] = MASK;

                int NNEig = GetNNeigbours(ip, NeigArr.data(), ndim, sizes.data());

                for(int i=0; i<NNEig; i++) {
                    if(Label[NeigArr[i]] < 0 && Label[NeigArr[i]] != WSHED)
                        continue;

                    Dist[ip] = 1;  //Initialise queue with neighbours at level h of current basins or watersheds
                    fifoQueue.push(ip);
                    break;
                }
            }

            int curdist = 1;
            fifoQueue.push(FICTITIOUS);

            while(true) { // extend basins
                int voxelIndex = fifoQueue.front();
                fifoQueue.pop();

                if(voxelIndex == FICTITIOUS) {
                    if(fifoQueue.empty())
                        break;

                    fifoQueue.push(FICTITIOUS);
                    curdist++;
                    voxelIndex = fifoQueue.front();
                    fifoQueue.pop();
                }

                int NNEig = GetNNeigbours(voxelIndex, NeigArr.data(), ndim, sizes.data());
                for(int i=0; i<NNEig; i++) { // Labelling p by inspecting neighbours
                    if(Dist[NeigArr[i]] < curdist && (Label[NeigArr[i]] > 0 || Label[NeigArr[i]]==WSHED)) {
                        if(Label[NeigArr[i]] > 0) { // q belongs to an existing basin or to a watershed
                            if(Label[voxelIndex] == MASK || Label[voxelIndex] ==WSHED)
                                Label[voxelIndex] = Label[NeigArr[i]]; // Removed from original algorithm || p.isLabelWSHED() )
                            else if(Label[voxelIndex] != Label[NeigArr[i]])
                                Label[voxelIndex] = WSHED;

                        } // end if lab>0
                        else if (Label[voxelIndex]==MASK)
                            Label[voxelIndex] = WSHED;
                    } else if(Label[NeigArr[i]]==MASK && Dist[NeigArr[i]]==0) {
                        Dist[NeigArr[i]] = curdist + 1;   //q is plateau pixel
                        fifoQueue.push(NeigArr[i]);
                    }
                } // end for, end processing neighbours
            } // end while (loop)

            // Detect and process new minima at level h
            for(int pixelIndex = CHist[h]; pixelIndex < CHist[h+1]; pixelIndex++) { //mask all pixels at level h
                int ip   = Index[pixelIndex];
                Dist[ip] = 0;       // Reset distance to zero

                if(Label[ip]!=MASK)
                    continue;
                curlab++;       // The pixel is inside a new minimum , create new label
                fifoQueue.push(ip);
                Label[ip] = curlab;

                while(fifoQueue.size()) {
                    int voxelIndex = fifoQueue.front();
                    fifoQueue.pop();

                    int NNEig = GetNNeigbours(voxelIndex, NeigArr.data(), ndim, sizes.data());  // replaced ip by voxelIndex
                    for(int i=0; i<NNEig; i++) { // inspect neighbours of q
                        if(Label[NeigArr[i]]!=MASK)
                            continue;

                        fifoQueue.push(NeigArr[i]);
                        Label[NeigArr[i]] = curlab;
                    }
                } // end while
            } // end for
        } // loop over h

        int MINS = (1<<15) -1;
        for ( unsigned i=0; i<NP; i++)
            data[i] = short(fMinOf2<int>(MINS, Label[i]));

        return true;
    }

}; // class

/** \brief overloaded operators +, *, - and / for basic numerical types
  *
  * This makes it possible to not only do imageNd +  <type>
  * but also                               <type> + imageNd, etc.
  */
template <typename T, typename U>
inline imageNd <T> operator+ ( U x, imageNd <T> y) {
    return y             + x;
}
template <typename T, typename U>
inline imageNd <T> operator* ( U x, imageNd <T> y) {
    return y             * x;
}
template <typename T, typename U>
inline imageNd <T> operator- ( U x, imageNd <T> y) {
    return -y            + x;
}
template <typename T, typename U>
inline imageNd <T> operator/ ( U x, imageNd <T> y) {
    return reciprocal(y) * x;
}

}; // namespace

#endif // MDIMAGE_HPP_INCLUDED
