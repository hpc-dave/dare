/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef UTILITIES_PROPERTYINFORMATION_H_
#define UTILITIES_PROPERTYINFORMATION_H_

namespace dare {

/*!
 * @brief a structure for communicating information at compile time
 * @tparam TagType type of the tag (usually an enum class)
 * @tparam Flag true or false
 * @tparam Tag the actual tag, part of TagType
 *
 * Consider following example:
 * enum class MyTag{ foo, bar, baz};
 * using MyInfo = FlaggedInfo<true, MyTag, MyTag::foo>;
 * using MyInfo2 = FlaggedInfo<false, MyTag, MyTag::bar>;
 * if constexpr(MyInfo::flag) {
 *   DoA();
 * }
 * if constexpr(!MyInfo2::flag) {
 *   DoB();
 * }
 */
template <bool Flag, typename TagType, TagType Tag>
struct FlaggedInfo {
    static const bool flag = Flag;
    using tag_type = TagType;
    static const tag_type tag = Tag;
};

/*! \struct TaggedTypeInfo
 * @tparam Type the type that should be communicated
 * @tparam TagType the type of the tag (usually an enum class)
 * @tparam Tag the tag, part of TagType
 * @brief a structure for communicating information about types at compile time
 * 
 * Consider following example:
 * enum class MyTag{ foo, bar, baz};
 * using MyInfo = TaggedTypeInfo<int, MyTag, MyTag::foo>;
 * using MyInfo2 = TaggedTypeInfo<double, MyTag, MyTag::bar>;
 * if constexpr(std::is_same_v<typename MyInfo::type, int>) {
 *  DoA();
 * }
 * static_assert(std::is_same_v<typename MyInfo2::type, void>, "invalid type provided");
 */
template<typename Type, typename TagType, TagType Tag>
struct TaggedTypeInfo {
    using type = Type;
    using tag_type = TagType;
    static const tag_type tag = Tag;
};

/*! \struct TaggedCountedTypeInfo
 * @tparam Type the type that should be communicated
 * @tparam TagType the type of the tag (usually an enum class)
 * @tparam Tag the tag, part of TagType
 * @tparam NUM_ENTITIES the number of entities of this type
 * @brief a structure for communicating information about types with a fixed number of entities at compile time
 */
template<typename Type, typename TagType, TagType Tag, std::size_t NUM_ENTITIES>
struct TaggedCountedTypeInfo : TaggedTypeInfo<Type, TagType, Tag> {
    static const std::size_t N = NUM_ENTITIES;
};

/*!
 * \brief a concept for flagged information types
 */
template <typename T>
concept FlaggedInfoType =
    requires {
        typename T::type;
        typename T::tag_type;
    } &&
    std::same_as<std::remove_cv_t<decltype(T::flag)>, bool>;

template <typename T>
struct is_flagged_info : std::false_type {};

template <FlaggedInfoType T>
struct is_flagged_info<T> : std::true_type {};

template <typename T>
constexpr bool is_flagged_info_v = is_flagged_info<T>::value;

/*!
 * \brief a concept for tagged information types
 */
template <typename T>
concept TaggedTypeInfoType =
    requires {
        typename T::type;
        typename T::tag_type;
    } &&
    std::same_as<std::remove_cv_t<decltype(T::tag)>, typename T::tag_type>;

template <typename T>
struct is_tagged_type_info : std::false_type {};

template <TaggedTypeInfoType T>
struct is_tagged_type_info<T> : std::true_type {};

template <typename T>
constexpr bool is_tagged_type_info_v = is_tagged_type_info<T>::value;

/*!
 * \brief a concept for tagged information types with a fixed number of entities
 */
template <typename T>
concept TaggedCountedTypeInfoType =
    TaggedTypeInfoType<T> &&
    std::same_as<std::remove_cv_t<typename T::N>, std::size_t>;

template <typename T>
struct is_tagged_counted_type_info : std::false_type {};

template <TaggedCountedTypeInfoType T>
struct is_tagged_counted_type_info<T> : std::true_type {};

template <typename T>
constexpr bool is_tagged_counted_type_info_v = is_tagged_counted_type_info<T>::value;

template<typename T, typename TagType, TagType Tag>
struct default_convert_to_tagged_info {
    using type = TaggedTypeInfo<T, TagType, Tag>;
};

template <TaggedTypeInfoType T, typename TagType, TagType Tag>
struct default_convert_to_tagged_info<T, TagType, Tag> {
    using type = T;
};

template <typename T, typename TagType, TagType Tag>
using default_convert_to_tagged_info_t = typename default_convert_to_tagged_info<T, TagType, Tag>::type;


/*!
 * \brief an tagging class for a None property
 */
struct None {
};

/*!
 * \brief concept to test for None types
 *
 * @tparam T type to test
 */
template <typename T>
concept NoneType =
    std::is_base_of_v<None, std::remove_cv_t<T>>;

/*!
 * \brief SFINAE test
 *
 * @tparam T type to test
 */
template <typename T>
struct is_none : std::bool_constant<NoneType<T>> {
};

template <typename T>
constexpr bool is_none_v = is_none<T>::value;

/*!
 * \brief an tagging class for a Multiple property
 */
struct Multiple {
};

/*!
 * \brief concept to test for a Multiple type
 *
 * @tparam T type to test
 */
template <typename T>
concept MultipleType =
    std::is_base_of_v<Multiple, std::remove_cv_t<T>>;

/*!
 * \brief SFINAE test
 *
 * @tparam T type to test
 */
template <typename T>
struct is_multiple : std::bool_constant<MultipleType<T>> {
};

template <typename T>
constexpr bool is_multiple_v = is_multiple<T>::value;

/*!
 * \brief an tagging class for a ConstCount property
 */
struct ConstCount {
};

/*!
 * \brief concept to test for a type
 *
 * @tparam T type to test
 */
template <typename T>
concept ConstCountType =
    std::is_base_of_v<ConstCount, std::remove_cv_t<T>>;

/*!
 * \brief SFINAE test
 *
 * @tparam T type to test
 */
template <typename T>
struct is_const_count : std::bool_constant<ConstCountType<T>> {
};

template <typename T>
constexpr bool is_const_count_v = is_const_count<T>::value;

template <typename T>
concept IntegerNumber = std::integral<decltype(T::value)> && std::is_convertible_v<T, int>;

template <typename T>
concept NaturalNumber =
    IntegerNumber<T>
    && requires {
        T::value >= 0;
};

template <std::size_t N>
using StaticNumber_size_t = std::integral_constant<std::size_t, N>;

static const StaticNumber_size_t<0> ZERO;
static const StaticNumber_size_t<1> ONE;
static const StaticNumber_size_t<2> TWO;
static const StaticNumber_size_t<3> THREE;
static const StaticNumber_size_t<4> FOUR;
static const StaticNumber_size_t<5> FIVE;
static const StaticNumber_size_t<6> SIX;
static const StaticNumber_size_t<7> SEVEN;
static const StaticNumber_size_t<8> EIGHT;
static const StaticNumber_size_t<9> NINE;
static const StaticNumber_size_t<10> TEN;

}  // namespace dare


// lets add some static tests
static_assert(dare::ZERO < dare::ONE);
static_assert(dare::ONE <= dare::ONE);
static_assert(dare::TEN > dare::TWO);
static_assert(dare::TEN >= dare::TWO);
static_assert(dare::THREE == 3);
static_assert(dare::IntegerNumber<decltype(dare::ZERO)>);
static_assert(dare::NaturalNumber<decltype(dare::ONE)>);
#endif  // UTILITIES_PROPERTYINFORMATION_H_
