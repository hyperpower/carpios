// This software is in the public domain. Where that dedication is not
// recognized, you are granted a perpetual, irrevocable license to copy,
// distribute, and modify this file as you see fit.
// Authored in 2015 by Dimitri Diakopoulos (http://www.dimitridiakopoulos.com)
// https://github.com/ddiakopoulos/tinyply

#ifndef _TS_IO_PLY_H_
#define _TS_IO_PLY_H_

#include <vector>
#include <algorithm>
#include <string>
#include <stdint.h>
#include <map>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <memory>
#include <functional>
#include <cstring>

namespace TS {

template<typename T>
T endian_swap(const T & v) {
	return v;
}
template<>
inline uint16_t endian_swap(const uint16_t & v) {
	return (v << 8) | (v >> 8);
}
template<>
inline uint32_t endian_swap(const uint32_t & v) {
	return (v << 24) | ((v << 8) & 0x00ff0000) | ((v >> 8) & 0x0000ff00)
			| (v >> 24);
}
template<>
inline uint64_t endian_swap(const uint64_t & v) {
	return (((v & 0x00000000000000ffLL) << 56)
			| ((v & 0x000000000000ff00LL) << 40)
			| ((v & 0x0000000000ff0000LL) << 24)
			| ((v & 0x00000000ff000000LL) << 8)
			| ((v & 0x000000ff00000000LL) >> 8)
			| ((v & 0x0000ff0000000000LL) >> 24)
			| ((v & 0x00ff000000000000LL) >> 40)
			| ((v & 0xff00000000000000LL) >> 56));
}
template<>
inline int16_t endian_swap(const int16_t & v) {
	uint16_t r = endian_swap(*(uint16_t*) &v);
	return *(int16_t*) &r;
}
template<>
inline int32_t endian_swap(const int32_t & v) {
	uint32_t r = endian_swap(*(uint32_t*) &v);
	return *(int32_t*) &r;
}
template<>
inline int64_t endian_swap(const int64_t & v) {
	uint64_t r = endian_swap(*(uint64_t*) &v);
	return *(int64_t*) &r;
}
inline float endian_swap_float(const uint32_t & v) {
	uint32_t r = endian_swap(v);
	return *(float*) &r;
}
inline double endian_swap_double(const uint64_t & v) {
	uint64_t r = endian_swap(v);
	return *(double*) &r;
}

struct DataCursor {
	void * vector;
	uint8_t * data;
	size_t offset;
	bool realloc = false;
};

class PlyProperty {
public:
	enum class Type
		: uint8_t
		{
			INVALID, INT8, UINT8, INT16, UINT16, INT32, UINT32, FLOAT32, FLOAT64
	};

	Type listType, propertyType;
	bool isList;
	int listCount = 0;
	std::string name;

	PlyProperty(std::istream & is) :
			isList(false) {
		parse_internal(is);
	}
	PlyProperty(Type type, const std::string & name) :
			propertyType(type), isList(false), name(name) {
	}
	PlyProperty(Type list_type, Type prop_type, const std::string & name,
			int listCount) :
			listType(list_type), propertyType(prop_type), isList(true), name(
					name), listCount(listCount) {
	}
protected:
	Type property_type_from_string(const std::string & t) {
		if (t == "int8" || t == "char")
			return Type::INT8;
		else if (t == "uint8" || t == "uchar")
			return Type::UINT8;
		else if (t == "int16" || t == "short")
			return Type::INT16;
		else if (t == "uint16" || t == "ushort")
			return Type::UINT16;
		else if (t == "int32" || t == "int")
			return Type::INT32;
		else if (t == "uint32" || t == "uint")
			return Type::UINT32;
		else if (t == "float32" || t == "float")
			return Type::FLOAT32;
		else if (t == "float64" || t == "double")
			return Type::FLOAT64;
		return Type::INVALID;
	}

	void parse_internal(std::istream & is) {
		std::string type;
		is >> type;
		if (type == "list") {
			std::string countType;
			is >> countType >> type;
			listType = property_type_from_string(countType);
			isList = true;
		}
		propertyType = property_type_from_string(type);
		is >> name;
	}

};

inline std::string make_key(const std::string & a, const std::string & b) {
	return (a + "-" + b);
}

template<typename T>
void ply_cast(void * dest, const char * src, bool be) {
	*(static_cast<T *>(dest)) =
			(be) ? endian_swap(*(reinterpret_cast<const T *>(src))) : *(reinterpret_cast<const T *>(src));
}

template<typename T>
void ply_cast_float(void * dest, const char * src, bool be) {
	*(static_cast<T *>(dest)) =
			(be) ? endian_swap_float(
							*(reinterpret_cast<const uint32_t *>(src))) :
					*(reinterpret_cast<const T *>(src));
}

template<typename T>
void ply_cast_double(void * dest, const char * src, bool be) {
	*(static_cast<T *>(dest)) =
			(be) ? endian_swap_double(
							*(reinterpret_cast<const uint64_t *>(src))) :
					*(reinterpret_cast<const T *>(src));
}

template<typename T>
T ply_read_ascii(std::istream & is) {
	T data;
	is >> data;
	return data;
}

template<typename T>
void ply_cast_ascii(void * dest, std::istream & is) {
	*(static_cast<T *>(dest)) = ply_read_ascii<T>(is);
}

struct PropertyInfo {
	int stride;
	std::string str;
};
static std::map<PlyProperty::Type, PropertyInfo> PropertyTable { {
		PlyProperty::Type::INT8, { 1, "char" } }, { PlyProperty::Type::UINT8, {
		1, "uchar" } }, { PlyProperty::Type::INT16, { 2, "short" } }, {
		PlyProperty::Type::UINT16, { 2, "ushort" } }, {
		PlyProperty::Type::INT32, { 4, "int" } }, { PlyProperty::Type::UINT32, {
		4, "uint" } }, { PlyProperty::Type::FLOAT32, { 4, "float" } }, {
		PlyProperty::Type::FLOAT64, { 8, "double" } }, {
		PlyProperty::Type::INVALID, { 0, "INVALID" } } };

template<typename T>
inline uint8_t * resize(void * v, size_t newSize) {
	auto vec = static_cast<std::vector<T> *>(v);
	vec->resize(newSize);
	return reinterpret_cast<uint8_t *>(vec->data());
}

inline void resize_vector(const PlyProperty::Type t, void * v, size_t newSize,
		uint8_t *& ptr) {
	switch (t) {
	case PlyProperty::Type::INT8:
		ptr = resize<int8_t>(v, newSize);
		break;
	case PlyProperty::Type::UINT8:
		ptr = resize<uint8_t>(v, newSize);
		break;
	case PlyProperty::Type::INT16:
		ptr = resize<int16_t>(v, newSize);
		break;
	case PlyProperty::Type::UINT16:
		ptr = resize<uint16_t>(v, newSize);
		break;
	case PlyProperty::Type::INT32:
		ptr = resize<int32_t>(v, newSize);
		break;
	case PlyProperty::Type::UINT32:
		ptr = resize<uint32_t>(v, newSize);
		break;
	case PlyProperty::Type::FLOAT32:
		ptr = resize<float>(v, newSize);
		break;
	case PlyProperty::Type::FLOAT64:
		ptr = resize<double>(v, newSize);
		break;
	case PlyProperty::Type::INVALID:
		throw std::invalid_argument("invalid ply property");
	}
}

template<typename T>
inline PlyProperty::Type property_type_for_type(std::vector<T> & theType) {
	if (std::is_same<T, int8_t>::value)
		return PlyProperty::Type::INT8;
	else if (std::is_same<T, uint8_t>::value)
		return PlyProperty::Type::UINT8;
	else if (std::is_same<T, int16_t>::value)
		return PlyProperty::Type::INT16;
	else if (std::is_same<T, uint16_t>::value)
		return PlyProperty::Type::UINT16;
	else if (std::is_same<T, int32_t>::value)
		return PlyProperty::Type::INT32;
	else if (std::is_same<T, uint32_t>::value)
		return PlyProperty::Type::UINT32;
	else if (std::is_same<T, float>::value)
		return PlyProperty::Type::FLOAT32;
	else if (std::is_same<T, double>::value)
		return PlyProperty::Type::FLOAT64;
	else
		return PlyProperty::Type::INVALID;
}

class PlyElement {
public:
	std::string name;
	size_t size;
	std::vector<PlyProperty> properties;

	PlyElement(std::istream & istream) {
		parse_internal(istream);
	}
	PlyElement(const std::string & name, size_t count) :
			name(name), size(count) {
	}

protected:
	void parse_internal(std::istream & is) {
		is >> name >> size;
	}
};

inline int find_element(const std::string key, std::vector<PlyElement> & list) {
	for (size_t i = 0; i < list.size(); ++i) {
		if (list[i].name == key) {
			return i;
		}
	}
	return -1;
}

class PlyFile {

public:

	bool isBinary = false;
	bool isBigEndian = false;

	std::map<std::string, std::shared_ptr<DataCursor> > userDataTable;

	std::vector<PlyElement> elements;
	std::vector<std::string> requestedElements;

	std::vector<std::string> comments;
	std::vector<std::string> objInfo;

	PlyFile() {

	}
	PlyFile(std::istream & is) {
		if (!parse_header(is)) {
			throw std::runtime_error(
					"file is not ply or encounted junk in header");
		}
	}

	void read(std::istream & is) {
		read_internal(is);
	}
	void write(std::ostream & os, bool isBinary) {
		if (isBinary) {
			write_binary_internal(os);
		} else {
			write_ascii_internal(os);
		}
	}

	std::vector<PlyElement> & get_elements() {
		return elements;
	}

	template<typename T>
	size_t request_properties_from_element( //
			const std::string & elementKey, //
			std::vector<std::string> propertyKeys, //
			std::vector<T>& source,               //
			const int listCount = 1) {             //
		if (elements.size() == 0)
			return 0;

		if (find_element(elementKey, elements) >= 0) {
			if (std::find(requestedElements.begin(), requestedElements.end(),
					elementKey) == requestedElements.end())
				requestedElements.push_back(elementKey);
		} else
			return 0;

		// count and verify large enough
		auto instance_counter =
				[&](const std::string & elementKey, const std::string & propertyKey)
				{
					for (auto e : elements)
					{
						if (e.name != elementKey) continue;
						for (auto p : e.properties)
						{
							if (p.name == propertyKey)
							{
								if (PropertyTable[property_type_for_type(source)].stride != PropertyTable[p.propertyType].stride)
								throw std::runtime_error("destination vector is wrongly typed to hold this property");
								return e.size;

							}
						}
					}
					return size_t(0);
				};

		// Check if requested key is in the parsed header
		std::vector<std::string> unusedKeys;
		for (auto key : propertyKeys) {
			for (auto e : elements) {
				if (e.name != elementKey)
					continue;
				std::vector<std::string> headerKeys;
				for (auto p : e.properties) {
					headerKeys.push_back(p.name);
				}

				if (std::find(headerKeys.begin(), headerKeys.end(), key)
						== headerKeys.end()) {
					unusedKeys.push_back(key);
				}

			}
		}

		// Not using them?
		// Don't let them affect the propertyKeys count used for calculating array sizes
		for (auto k : unusedKeys) {
			propertyKeys.erase(
					std::remove(propertyKeys.begin(), propertyKeys.end(), k),
					propertyKeys.end());
		}
		if (!propertyKeys.size())
			return 0;

		// All requested properties in the userDataTable share the same cursor (thrown into the same flat array)
		auto cursor = std::make_shared<DataCursor>();

		std::vector<size_t> instanceCounts;

		for (auto key : propertyKeys) {
			if (int instanceCount = instance_counter(elementKey, key)) {
				instanceCounts.push_back(instanceCount);
				auto result = userDataTable.insert(
						std::pair<std::string, std::shared_ptr<DataCursor>>(
								make_key(elementKey, key), cursor));
				if (result.second == false)
					throw std::invalid_argument(
							"property has already been requested: " + key);
			} else
				continue;
		}

		size_t totalInstanceSize = [&]() {
			size_t t = 0;
			for (auto c : instanceCounts) {
				t += c;
			}
			return t;}() * listCount;
		source.resize(totalInstanceSize); // this satisfies regular properties; `cursor->realloc` is for list types since tinyply uses single-pass parsing
		cursor->offset = 0;
		cursor->vector = &source;
		cursor->data = reinterpret_cast<uint8_t *>(source.data());

		if (listCount > 1) {
			cursor->realloc = true;
			return (totalInstanceSize / propertyKeys.size()) / listCount;
		}

		return totalInstanceSize / propertyKeys.size();
	}

	template<typename T>
	void add_properties_to_element(const std::string & elementKey,
			const std::vector<std::string> & propertyKeys,
			std::vector<T> & source, const int listCount = 1,
			const PlyProperty::Type listType = PlyProperty::Type::INVALID) {
		auto cursor = std::make_shared<DataCursor>();
		cursor->offset = 0;
		cursor->vector = &source;
		cursor->data = reinterpret_cast<uint8_t *>(source.data());

		auto create_property_on_element =
				[&](PlyElement & e)
				{
					for (auto key : propertyKeys)
					{
						PlyProperty::Type t = property_type_for_type(source);
						PlyProperty newProp =
						(listType == PlyProperty::Type::INVALID) ?
						PlyProperty(t, key) :
						PlyProperty(listType, t, key, listCount);
						userDataTable.insert(
								std::pair<std::string, std::shared_ptr<DataCursor> >(make_key(e.name, key), cursor));
						e.properties.push_back(newProp);
					}
				};

		int idx = find_element(elementKey, elements);
		if (idx >= 0) {
			PlyElement & e = elements[idx];
			create_property_on_element(e);
		} else {
			PlyElement newElement =
					(listCount == 1) ?
							PlyElement(elementKey,
									source.size() / propertyKeys.size()) :
							PlyElement(elementKey, source.size() / listCount);
			create_property_on_element(newElement);
			elements.push_back(newElement);
		}
	}

private:

	size_t skip_property_binary(const PlyProperty & property,
			std::istream & is) {
		static std::vector<char> skip(
				PropertyTable[property.propertyType].stride);
		if (property.isList) {
			size_t listSize = 0;
			size_t dummyCount = 0;
			read_property_binary(property.listType, &listSize, dummyCount, is);
			for (size_t i = 0; i < listSize; ++i)
				is.read(skip.data(),
						PropertyTable[property.propertyType].stride);
			return listSize;
		} else {
			is.read(skip.data(), PropertyTable[property.propertyType].stride);
			return 0;
		}
	}
	void skip_property_ascii(const PlyProperty & property, std::istream & is) {
		std::string skip;
		if (property.isList) {
			int listSize;
			is >> listSize;
			for (int i = 0; i < listSize; ++i) {
				is >> skip;
			}
		} else {
			is >> skip;
		}
	}

	void read_property_binary(PlyProperty::Type t, void * dest,
			size_t & destOffset, std::istream & is) {
		static std::vector<char> src(PropertyTable[t].stride);
		is.read(src.data(), PropertyTable[t].stride);

		switch (t) {
		case PlyProperty::Type::INT8:
			ply_cast<int8_t>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::UINT8:
			ply_cast<uint8_t>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::INT16:
			ply_cast<int16_t>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::UINT16:
			ply_cast<uint16_t>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::INT32:
			ply_cast<int32_t>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::UINT32:
			ply_cast<uint32_t>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::FLOAT32:
			ply_cast_float<float>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::FLOAT64:
			ply_cast_double<double>(dest, src.data(), isBigEndian);
			break;
		case PlyProperty::Type::INVALID:
			throw std::invalid_argument("invalid ply property");
		}
		destOffset += PropertyTable[t].stride;
	}

	void read_property_ascii(PlyProperty::Type t, void * dest,
			size_t& destOffset, std::istream & is) {
		switch (t) {
		case PlyProperty::Type::INT8:
			*((int8_t *) dest) = ply_read_ascii<int32_t>(is);
			break;
		case PlyProperty::Type::UINT8:
			*((uint8_t *) dest) = ply_read_ascii<uint32_t>(is);
			break;
		case PlyProperty::Type::INT16:
			ply_cast_ascii<int16_t>(dest, is);
			break;
		case PlyProperty::Type::UINT16:
			ply_cast_ascii<uint16_t>(dest, is);
			break;
		case PlyProperty::Type::INT32:
			ply_cast_ascii<int32_t>(dest, is);
			break;
		case PlyProperty::Type::UINT32:
			ply_cast_ascii<uint32_t>(dest, is);
			break;
		case PlyProperty::Type::FLOAT32:
			ply_cast_ascii<float>(dest, is);
			break;
		case PlyProperty::Type::FLOAT64:
			ply_cast_ascii<double>(dest, is);
			break;
		case PlyProperty::Type::INVALID:
			throw std::invalid_argument("invalid ply property");
		}
		destOffset += PropertyTable[t].stride;
	}
	void write_property_ascii(PlyProperty::Type t, std::ostream & os,
			uint8_t * src, size_t & srcOffset) {
		switch (t) {
		case PlyProperty::Type::INT8:
			os << static_cast<int32_t>(*reinterpret_cast<int8_t*>(src));
			break;
		case PlyProperty::Type::UINT8:
			os << static_cast<uint32_t>(*reinterpret_cast<uint8_t*>(src));
			break;
		case PlyProperty::Type::INT16:
			os << *reinterpret_cast<int16_t*>(src);
			break;
		case PlyProperty::Type::UINT16:
			os << *reinterpret_cast<uint16_t*>(src);
			break;
		case PlyProperty::Type::INT32:
			os << *reinterpret_cast<int32_t*>(src);
			break;
		case PlyProperty::Type::UINT32:
			os << *reinterpret_cast<uint32_t*>(src);
			break;
		case PlyProperty::Type::FLOAT32:
			os << *reinterpret_cast<float*>(src);
			break;
		case PlyProperty::Type::FLOAT64:
			os << *reinterpret_cast<double*>(src);
			break;
		case PlyProperty::Type::INVALID:
			throw std::invalid_argument("invalid ply property");
		}
		os << " ";
		srcOffset += PropertyTable[t].stride;
	}
	void write_property_binary(PlyProperty::Type t, std::ostream & os,
			uint8_t * src, size_t & srcOffset) {
		os.write((char *) src, PropertyTable[t].stride);
		srcOffset += PropertyTable[t].stride;
	}

	bool parse_header(std::istream & is) {
		std::string line;
		bool gotMagic = false;
		while (std::getline(is, line)) {
			std::istringstream ls(line);
			std::string token;
			ls >> token;
			if (token == "ply" || token == "PLY" || token == "") {
				gotMagic = true;
				continue;
			} else if (token == "comment")
				read_header_text(line, ls, comments, 8);
			else if (token == "format")
				read_header_format(ls);
			else if (token == "element")
				read_header_element(ls);
			else if (token == "property")
				read_header_property(ls);
			else if (token == "obj_info")
				read_header_text(line, ls, objInfo, 9);
			else if (token == "end_header")
				break;
			else
				return false;
		}
		return true;
	}
	void write_header(std::ostream & os) {
		const std::locale & fixLoc = std::locale("C");
		os.imbue(fixLoc);

		os << "ply" << std::endl;
		if (isBinary)
			os
					<< ((isBigEndian) ?
							"format binary_big_endian 1.0" :
							"format binary_little_endian 1.0") << std::endl;
		else
			os << "format ascii 1.0" << std::endl;

		for (const auto & comment : comments)
			os << "comment " << comment << std::endl;

		for (auto & e : elements) {
			os << "element " << e.name << " " << e.size << std::endl;
			for (const auto & p : e.properties) {
				if (p.isList) {
					os << "property list " << PropertyTable[p.listType].str
							<< " " << PropertyTable[p.propertyType].str << " "
							<< p.name << std::endl;
				} else {
					os << "property " << PropertyTable[p.propertyType].str
							<< " " << p.name << std::endl;
				}
			}
		}
		os << "end_header" << std::endl;
	}

	void read_header_format(std::istream & is) {
		std::string s;
		(is >> s);
		if (s == "binary_little_endian")
			isBinary = true;
		else if (s == "binary_big_endian")
			isBinary = isBigEndian = true;
	}
	void read_header_element(std::istream & is) {
		elements.emplace_back(is);
	}
	void read_header_property(std::istream & is) {
		elements.back().properties.emplace_back(is);
	}
	void read_header_text(std::string line, std::istream & is,
			std::vector<std::string> & place, int erase = 0) {
		place.push_back((erase > 0) ? line.erase(0, erase) : line);
	}

	void read_internal(std::istream & is) {
		std::function<
				void(PlyProperty::Type t, void * dest, size_t & destOffset,
						std::istream & is)> read;
		std::function<void(const PlyProperty & property, std::istream & is)> skip;
		if (isBinary) {
			read =
					[&](PlyProperty::Type t, void * dest, size_t & destOffset, std::istream & is) {
						read_property_binary(t, dest, destOffset, is);
					};
			skip =
					[&](const PlyProperty & property, std::istream & is) {skip_property_binary(property, is);};
		} else {
			read =
					[&](PlyProperty::Type t, void * dest, size_t & destOffset, std::istream & is) {
						read_property_ascii(t, dest, destOffset, is);
					};
			skip =
					[&](const PlyProperty & property, std::istream & is) {skip_property_ascii(property, is);};
		}

		for (auto & element : get_elements()) {
			if (std::find(requestedElements.begin(), requestedElements.end(),
					element.name) != requestedElements.end()) {
				//std::cout <<"name : "<< element.name << "\n";
				for (size_t count = 0; count < element.size; ++count) {
					for (auto & property : element.properties) {
						//std::cout <<"property : "<< property.name << "\n";
						auto & cursor = userDataTable[make_key(element.name,
								property.name)];
						if (cursor != nullptr) {
							if (property.isList) {
								size_t listSize = 0;
								size_t dummyCount = 0;
								read(property.listType,
										&listSize,
										dummyCount,
										is);
								if (cursor->realloc == false) {
									cursor->realloc = true;
									resize_vector(property.propertyType,
											cursor->vector,
											listSize * element.size,
											cursor->data);
								}
								for (size_t i = 0; i < listSize; ++i) {
									read(property.propertyType,
										(cursor->data + cursor->offset),
										cursor->offset, is);
								}
							} else {
								read(property.propertyType,
										(cursor->data + cursor->offset),
										cursor->offset, is);
							}
						} else {
							skip(property, is);
						}
					}
				}
			} else
				continue;
		}
	}

	void write_ascii_internal(std::ostream & os) {
		write_header(os);

		for (auto & e : elements) {
			for (size_t i = 0; i < e.size; ++i) {
				for (auto & p : e.properties) {
					auto & cursor = userDataTable[make_key(e.name, p.name)];
					if (p.isList) {
						os << p.listCount << " ";
						for (int j = 0; j < p.listCount; ++j) {
							write_property_ascii(p.propertyType, os,
									(cursor->data + cursor->offset),
									cursor->offset);
						}
					} else {
						write_property_ascii(p.propertyType, os,
								(cursor->data + cursor->offset),
								cursor->offset);
					}
				}
				os << std::endl;
			}
		}
	}

	void write_binary_internal(std::ostream & os) {
		isBinary = true;
		write_header(os);

		for (auto & e : elements) {
			for (size_t i = 0; i < e.size; ++i) {
				for (auto & p : e.properties) {
					auto & cursor = userDataTable[make_key(e.name, p.name)];
					if (p.isList) {
						uint8_t listSize[4] = { 0, 0, 0, 0 };
						memcpy(listSize, &p.listCount, sizeof(uint32_t));
						size_t dummyCount = 0;
						write_property_binary(p.listType, os, listSize,
								dummyCount);
						for (int j = 0; j < p.listCount; ++j) {
							write_property_binary(p.propertyType, os,
									(cursor->data + cursor->offset),
									cursor->offset);
						}
					} else {
						write_property_binary(p.propertyType, os,
								(cursor->data + cursor->offset),
								cursor->offset);
					}
				}
			}
		}
	}

};

} // namesapce ts

#endif // ts_io_ply_h
