#pragma once
#ifndef ANDRES_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_HDF5_HXX

#include <cassert>
#include <string>
#include <vector>
#include <array>
#include <sstream>

// compat fix for buggy hdf5 1.8 versions
#include <H5version.h>
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8 && defined(H5_USE_16_API_DEFAULT))
#define H5Gcreate_vers 2
#define H5Gopen_vers 2
#define H5Dopen_vers 2
#define H5Dcreate_vers 2
#define H5Acreate_vers 2
#endif

#include "hdf5.h"

#ifndef NDEBUG
#define ANDRES_GRAPH_HDF5_DEBUG true
#else
#define ANDRES_GRAPH_HDF5_DEBUG false
#endif

namespace andres{
namespace graph{
namespace hdf5{

template<class GRAPH>
class GraphTraitsHDF5;

// \cond suppress doxygen
template<bool B = true> class HandleCheck;
template<> class HandleCheck<true> {
public:
    HandleCheck() 
        { counter_ = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL); }
    ~HandleCheck()
        { assert( counter_ == H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL)); }
private:
    ssize_t counter_;
};

template<> class HandleCheck<false> {
public:
    void check() {}
};
// \endcond

template<class T>
inline hid_t uintTypeHelper() {
    switch(sizeof(T)) {
        case 1: return H5T_STD_U8LE;
        case 2: return H5T_STD_U16LE;
        case 4: return H5T_STD_U32LE;
        case 8: return H5T_STD_U64LE;
        default: throw std::runtime_error("No matching HDF5 type.");
    }
}

template<class T>
inline hid_t intTypeHelper() {
    switch(sizeof(T)) {
        case 1: return H5T_STD_I8LE;
        case 2: return H5T_STD_I16LE;
        case 4: return H5T_STD_I32LE;
        case 8: return H5T_STD_I64LE;
        default: throw std::runtime_error("No matching HDF5 type.");
    }
}

template<class T>
inline hid_t floatingTypeHelper() {
    switch(sizeof(T)) {
        case 4: return H5T_IEEE_F32LE;
        case 8: return H5T_IEEE_F64LE;
        default: throw std::runtime_error("No matching HDF5 type.");
    }
}

template<class T>
inline hid_t hdf5Type();

template<> inline hid_t hdf5Type<unsigned char>()
    { return uintTypeHelper<unsigned char>(); }
template<> inline hid_t hdf5Type<unsigned short>()
    { return uintTypeHelper<unsigned short>(); }
template<> inline hid_t hdf5Type<unsigned int>()
    { return uintTypeHelper<unsigned int>(); }
template<> inline hid_t hdf5Type<unsigned long>()
    { return uintTypeHelper<unsigned long>(); }
template<> inline hid_t hdf5Type<unsigned long long>()
    { return uintTypeHelper<unsigned long long>(); }
template<> inline hid_t hdf5Type<char>()
    { return uintTypeHelper<char>(); }

template<> inline hid_t hdf5Type<signed char>()
    { return intTypeHelper<signed char>(); }
template<> inline hid_t hdf5Type<short>()
    { return intTypeHelper<short>(); }
template<> inline hid_t hdf5Type<int>()
    { return intTypeHelper<int>(); }
template<> inline hid_t hdf5Type<long>()
    { return intTypeHelper<long>(); }
template<> inline hid_t hdf5Type<long long>()
    { return intTypeHelper<long long>(); }

template<> inline hid_t hdf5Type<float>()
    { return floatingTypeHelper<float>(); }
template<> inline hid_t hdf5Type<double>()
    { return floatingTypeHelper<double>(); }

enum class FileAccessMode {READ_ONLY, READ_WRITE};
enum class HDF5Version {DEFAULT_HDF5_VERSION, LATEST_HDF5_VERSION};

hid_t createFile(const std::string&, HDF5Version = HDF5Version::DEFAULT_HDF5_VERSION);
hid_t openFile(const std::string&, FileAccessMode = FileAccessMode::READ_ONLY, HDF5Version = HDF5Version::DEFAULT_HDF5_VERSION);
void closeFile(const hid_t&);

hid_t createGroup(const hid_t&, const std::string& groupName);
hid_t openGroup(const hid_t&, const std::string&, const bool =false);
void closeGroup(const hid_t&);

template<class T>
void save(const hid_t&, const std::string&, std::initializer_list<std::size_t>, const T * const);
template<class T>
void save(const hid_t&, const std::string&, const T&);

template<class T>
void load(const hid_t&, const std::string&, std::vector<std::size_t>&, std::vector<T>&);
template<class T>
void load(const hid_t&, const std::string&, std::vector<std::vector<T> >&);
template<class T>
void load(const hid_t&, const std::string&, T&);

/// Create an HDF5 file.
///
/// \param filename Name of the file.
/// \param hdf5version HDF5 version tag.
///
/// \returns HDF5 handle
///
/// \sa openFile(), closeFile()
///
inline hid_t
createFile (
    const std::string& filename,
    HDF5Version hdf5version
) {
    hid_t version = H5P_DEFAULT;
    if(hdf5version == HDF5Version::LATEST_HDF5_VERSION) {
        version = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_libver_bounds(version, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    }
    hid_t fileHandle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, version);
    H5Pclose(version);
    if(fileHandle < 0) {
        throw std::runtime_error("Could not create HDF5 file: " + filename);
    }
    return fileHandle;
}

/// Open an HDF5 file.
///
/// \param filename Name of the file.
/// \param fileAccessMode File access mode.
/// \param hdf5version HDF5 version tag.
///
/// \returns HDF5 handle
///
/// \sa closeFile(), createFile()
///
inline hid_t
openFile(
    const std::string& filename,
    FileAccessMode fileAccessMode,
    HDF5Version hdf5version
) {
    hid_t access = H5F_ACC_RDONLY;
    if(fileAccessMode == FileAccessMode::READ_WRITE) {
        access = H5F_ACC_RDWR;
    }

    hid_t version = H5P_DEFAULT;
    if(hdf5version == HDF5Version::LATEST_HDF5_VERSION) {
        version = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_libver_bounds(version, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    }

    hid_t fileHandle = H5Fopen(filename.c_str(), access, version);
    if(hdf5version == HDF5Version::LATEST_HDF5_VERSION) {
        H5Pclose(version);
    }
    if(fileHandle < 0) {
        throw std::runtime_error("Could not open HDF5 file: " + filename);
    }
    return fileHandle;
}

/// Close an HDF5 file
/// 
/// \param handle Handle to the HDF5 file.
///
/// \sa openFile(), createFile()
///
inline void closeFile (
    const hid_t& handle
) {
    H5Fclose(handle);
}

/// Create an HDF5 group.
///
/// \param parentHandle HDF5 handle on the parent group or file.
/// \param groupName Name of the group.
/// \returns HDF5 handle on the created group
///
/// \sa openGroup(), closeGroup()
///
inline hid_t 
createGroup (
    const hid_t& parentHandle,
    const std::string& groupName
) { 
    hid_t groupHandle = H5Gcreate(parentHandle, groupName.c_str(), 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(groupHandle < 0) {
        throw std::runtime_error("Could not create HDF5 group.");
    }
    return groupHandle;
}

/// Open an HDF5 group.
///
/// \param parentHandle HDF5 handle on the parent group or file.
/// \param groupName Name of the group.
/// \param createIfNonexistent if set to true, the group is newly created if it cannot be opened.
/// \returns HDF5 handle on the opened group.
///
/// \sa createGroup(), closeGroup()
///
inline hid_t 
openGroup (
    const hid_t& parentHandle,
    const std::string& groupName,
    const bool createIfNonexistent
) { 
    hid_t groupHandle;

    if(createIfNonexistent) {
        H5E_auto2_t old_func;
        void *old_client_data;
        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);
        groupHandle = H5Gopen(parentHandle, groupName.c_str(), H5P_DEFAULT);
        H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
        if(groupHandle != 0) {
            groupHandle = createGroup(parentHandle, groupName);
        }
        if(groupHandle < 0) {
            throw std::runtime_error("HDF5: Could not open or create HDF5 group " + groupName + ".");
        }
    } else {
        groupHandle = H5Gopen(parentHandle, groupName.c_str(), H5P_DEFAULT);
        if(groupHandle < 0) {
            throw std::runtime_error("HDF5: Could not open HDF5 group " + groupName + ".");
        }
    }
    return groupHandle;
}

/// Close an HDF5 group.
///
/// \param handle HDF5 handle on group to close.
///
/// \sa openGroup(), createGroup()
///
inline void 
closeGroup(
    const hid_t& handle
) {
    H5Gclose(handle);
}

/// Save shaped data as HDF5 dataset.
///
/// \param parentHandle Handle of the parent HDF5 file or group.
/// \param datasetName Name of the HDF5 dataset.
/// \param shape Shape of the data.
/// \param data A total of prod(shape) values must be readable.
///
template<class T>
void save(
    const hid_t& parentHandle,
    const std::string& datasetName,
    const std::initializer_list<std::size_t> shape,
    const T * const data
) {
    assert(parentHandle >= 0);
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;

    hid_t datatype, dataspace, dataset;
    std::string sError;
    
    // build dataspace
    datatype = H5Tcopy(hdf5Type<T>());
    {
        std::vector<hsize_t> storeShape(shape.begin(), shape.end());
        dataspace = H5Screate_simple(storeShape.size(), &storeShape[0], NULL);
    }
    if(dataspace < 0) {
        sError = "Dataspace creation failed.";
        goto cleanupType;
    }

    // create new dataset
    dataset = H5Dcreate(parentHandle, datasetName.c_str(), datatype,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dataset < 0) {
        sError = "Dataset creation failed.";
        goto cleanupSpace;
    }

    if(H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL,H5P_DEFAULT, data)<0) {
        sError = "Failed  to write data to dataset.";
        goto cleanupDataset;
    }

cleanupDataset:
    H5Dclose(dataset);
cleanupSpace:
    H5Sclose(dataspace);
cleanupType:
    H5Tclose(datatype);
    if(!sError.empty()) {
        throw std::runtime_error("HDF5: Saving sataset '"+datasetName+"' failed: " + sError);
    }
}

/// Save a scalar as an HDF5 dataset.
///
/// \param parentHandle Handle of the parent HDF5 file or group.
/// \param datasetName Name of the HDF5 dataset.
/// \param data Scalar to save.
///
template<class T>
void save(
    const hid_t& parentHandle,
    const std::string& datasetName,
    const T& data
) {
    assert(parentHandle >= 0);
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;

    std::string sError;
    hid_t datatype, dataspace, dataset;
    // build dataspace
    datatype = H5Tcopy(hdf5Type<T>());
    dataspace = H5Screate(H5S_SCALAR);
    if(dataspace < 0) {
        sError = "Cannot create dataspace.";
        goto cleanupType;
    }

    // create new dataset
    dataset = H5Dcreate(parentHandle, datasetName.c_str(), datatype,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dataset < 0) {
        sError = "Cannot create dataset.";
        goto cleanupSpace;
    }

    if(H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL,H5P_DEFAULT, &data)!=0) {
        sError = "Cannot write data.";
    }
    
    H5Dclose(dataset);

cleanupSpace:
    H5Sclose(dataspace);
cleanupType:
    H5Tclose(datatype);
    if(!sError.empty()) {
        throw std::runtime_error(
            "HDF5: Saving dataset '"+datasetName+"' failed: " + sError
        );
    }
}

/// Load a 2D array from an HDF5 dataset as a single vector.
///
/// \param parentHandle Handle of the parent HDF5 file or group.
/// \param datasetName Name of the HDF5 dataset.
/// \param shape The shape of the read data.
/// If an error occures, shape is cleared and an exception is thrown.
/// \param data Read data.
///
template<class T>
void load(
    const hid_t& parentHandle,
    const std::string& datasetName,
    std::vector<std::size_t>& shape,
    std::vector<T>& data
) {
    assert(parentHandle >= 0);
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;

    hid_t dataset, memspace, filespace, type, nativeType;
    int ndims;
    std::vector<hsize_t> storedShape;
    std::string sError;
    shape.clear();
    
    dataset = H5Dopen(parentHandle, datasetName.c_str(), H5P_DEFAULT);
    if(dataset < 0) {
        sError = "Cannot open dataset.";
        goto errorCheck;
    }
    filespace = H5Dget_space(dataset);
    type = H5Dget_type(dataset);
    nativeType = H5Tget_native_type(type, H5T_DIR_DESCEND);
    if(!H5Tequal(nativeType, hdf5Type<T>())) {
        sError = "Requested data type does not matched stored one.";
        goto cleanupDataset;
    }
    
    // Retrieve dimensions
    ndims = H5Sget_simple_extent_ndims(filespace);
    if(ndims<=0) {
        sError = ndims == 0?"Dataset is a scalar.":"Cannot get number of dimensions.";
        goto cleanupDataset;
    }
    storedShape.resize(ndims);
    if(H5Sget_simple_extent_dims(filespace, &storedShape[0], NULL)<0) {
        sError = "Cannot get simple extent dimensions.";
        goto cleanupDataset;
    }
    memspace = H5Screate_simple(ndims, &storedShape[0], NULL);
    if(memspace<=0) {
        sError = "Could not allocate memory space.";
        goto cleanupTypes;
    }
    
    // Read data
    {
        std::size_t numberOfElements = 1;
        for(std::size_t i=0;i<storedShape.size();++i)
            numberOfElements *= storedShape[i];
        data.resize(numberOfElements);
    }
    if(H5Dread(dataset, nativeType, memspace, filespace,H5P_DEFAULT, &data[0])!=0) {
        sError = "Could not read data from dataset.";
        goto cleanupDataset;
    }
    shape.assign(storedShape.begin(),storedShape.end());

cleanupDataset:
    H5Dclose(dataset);
    H5Sclose(memspace);
cleanupTypes:
    H5Tclose(nativeType);
    H5Tclose(type);
    H5Sclose(filespace);
errorCheck:
    if(!sError.empty()) {
        throw std::runtime_error("HDF5: Loading dataset '"+datasetName+"' failed: "+sError);
    }
}

/// Load a scalar value from an HDF5 dataset.
///
/// \param parentHandle Handle of the parent HDF5 file or group.
/// \param datasetName Name of the HDF5 dataset.
/// \param data Scalar.
///
/// \sa loadHyperslab()
///
template<class T>
void load(
    const hid_t& parentHandle,
    const std::string& datasetName,
    T& data
) {
    assert(parentHandle >= 0);
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    
    std::string sError;
    hid_t dataset, filespace, type, nativeType, memspace;
    int ndims;
    dataset = H5Dopen(parentHandle, datasetName.c_str(), H5P_DEFAULT);
    if(dataset < 0) {
        sError = "Cannot open dataset.";
        goto error;
    }
    filespace = H5Dget_space(dataset);
    type = H5Dget_type(dataset);
    nativeType = H5Tget_native_type(type, H5T_DIR_DESCEND);
    if(!H5Tequal(nativeType, hdf5Type<T>())) {
        sError = "Stored data type does not match the one requested.";
        goto cleanupTypes;
    }
    ndims = H5Sget_simple_extent_ndims(filespace);
    if(ndims!=0) {
        sError = "Dataset dimension is not zero.";
        goto cleanupTypes;
    }
    memspace = H5Screate(H5S_SCALAR);
    if(memspace<0) {
        sError = "Failed to create memory dataspace.";
        goto cleanupAll;
    }
    if(H5Dread(dataset, nativeType, memspace, filespace,H5P_DEFAULT, &data)<0) {
        sError = "Failed to write scalar.";
        goto cleanupAll;
    }
    
cleanupAll:
    H5Sclose(memspace);
cleanupTypes:
    H5Tclose(nativeType);
    H5Tclose(type);
    H5Sclose(filespace);
    H5Dclose(dataset);
error:
    if(!sError.empty()) {
        throw std::runtime_error(
            "HDF5:Loading scalar from dataset '"+datasetName+"' failed: " + sError
        );
    }
}


/// Load a string value from an HDF5 dataset.
///
/// \param parentHandle Handle of the parent HDF5 file or group.
/// \param datasetName Name of the HDF5 dataset.
/// \param data string to hold the read data.
///
/// \sa loadHyperslab()
///
template<>
inline void
load(
    const hid_t& parentHandle,
    const std::string& datasetName,
    std::string& data
) {
    assert(parentHandle>0);
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    
    std::string sError;
    
    hid_t dataSet, fileType, dataSpace, memType;
    hsize_t stringSize;
    
    dataSet = H5Dopen (parentHandle, datasetName.c_str(), H5P_DEFAULT);
    fileType = H5Dget_type(dataSet);
    dataSpace = H5Dget_space(dataSet);
    // Verify we opened an 1 dimensional array.
    {
        hsize_t dims[1];
        hsize_t ndims = H5Sget_simple_extent_dims(dataSpace, NULL, NULL);
        if(ndims!=1) {
            sError = "Dataset spans more than 1 dimensions.";
            goto cleanupTypes;
        }
        ndims = H5Sget_simple_extent_dims(dataSpace, dims, NULL);
        // Assume only 1 string saved
        if(dims[0]!=1) {
            sError = "Dataset contains more than 1 elements.";
            goto cleanupTypes;
        }
    }
    // Verify this is not a variable dimension array
    {
        htri_t isVariable = H5Tis_variable_str(fileType);
        if(isVariable) {
            sError = "String type is unsupported.";
            goto cleanupTypes;
            //const hid_t memType = H5Tget_native_type(fileType, H5T_DIR_ASCEND);
            //const hid_t mc = H5Tcopy(memType);
            // something like: char **ptrBuff = (char **) malloc(dim[0]*sizeof(char *));
            //status = H5Dread (dataSet, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ptrBuf);
            // Use
            //status = H5Dvlen_reclaim (memType, dataSpace, H5P_DEFAULT, ptrBuf);
        }
    }
    stringSize = H5Tget_size(fileType);
    memType = H5Tcopy(H5T_C_S1);
    H5Tset_size(memType,stringSize);
    {
        std::vector<char> buf(stringSize);
        herr_t status = H5Dread (dataSet, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0]);
        if(status!=0) {
            sError = "Could not read string data.";
            goto cleanupAll;
        }
        data.assign(&buf[0]);
    }

cleanupAll:    
    H5Tclose (memType);
cleanupTypes:
    H5Sclose (dataSpace);
    H5Tclose (fileType);
    H5Dclose (dataSet);
    if(!sError.empty()) {
        throw std::runtime_error(
            "HDF5: Loading string from dataset '"+datasetName+"' failed: "+sError
        );
    }
}


} //namespace hdf
} //namespace graph
} //namespace andres

#endif // #ifndef ANDRES_GRAPH_HDF5_HXX
