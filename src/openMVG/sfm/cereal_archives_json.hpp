#ifndef CEREAL_ARCHIVES_JSON_NO_PRETTY_WRITER_HPP_
#define CEREAL_ARCHIVES_JSON_NO_PRETTY_WRITER_HPP_

#include "cereal/cereal.hpp"
#include "cereal/details/util.hpp"

#include "cereal/external/rapidjson/writer.h"
#include "cereal/external/rapidjson/ostreamwrapper.h"
#include "cereal/external/rapidjson/istreamwrapper.h"
#include "cereal/external/rapidjson/document.h"
#include "cereal/external/base64.hpp"

#include <limits>
#include <sstream>
#include <stack>
#include <vector>
#include <string>

namespace cereal
{
  class JSONOutputArchiveNoPrettyWriter : public OutputArchive<JSONOutputArchiveNoPrettyWriter>, public traits::TextArchive
  {
    enum class NodeType { StartObject, InObject, StartArray, InArray };

    using WriteStream = CEREAL_RAPIDJSON_NAMESPACE::OStreamWrapper;
    using JSONWriter = CEREAL_RAPIDJSON_NAMESPACE::Writer<WriteStream>;

    public:
      /*! @name Common Functionality
          Common use cases for directly interacting with an JSONOutputArchiveNoPrettyWriter */
      //! @{

      //! A class containing various advanced options for the JSON archive
      class Options
      {
        public:
          //! Default options
          static Options Default(){ return Options(); }

          //! Default options with no indentation
          static Options NoIndent(){ return Options( JSONWriter::kDefaultMaxDecimalPlaces, IndentChar::space, 0 ); }

          //! The character to use for indenting
          enum class IndentChar : char
          {
            space = ' ',
            tab = '\t',
            newline = '\n',
            carriage_return = '\r'
          };

          //! Specify specific options for the JSONOutputArchiveNoPrettyWriter
          /*! @param precision The precision used for floating point numbers
              @param indentChar The type of character to indent with
              @param indentLength The number of indentChar to use for indentation
                             (0 corresponds to no indentation) */
          explicit Options( int precision = JSONWriter::kDefaultMaxDecimalPlaces,
                            IndentChar indentChar = IndentChar::space,
                            unsigned int indentLength = 4 ) :
            itsPrecision( precision ),
            itsIndentChar( static_cast<char>(indentChar) ),
            itsIndentLength( indentLength ) { }

        private:
          friend class JSONOutputArchiveNoPrettyWriter;
          int itsPrecision;
          char itsIndentChar;
          unsigned int itsIndentLength;
      };

      //! Construct, outputting to the provided stream
      /*! @param stream The stream to output to.
          @param options The JSON specific options to use.  See the Options struct
                         for the values of default parameters */
      JSONOutputArchiveNoPrettyWriter(std::ostream & stream, Options const & options = Options::Default() ) :
        OutputArchive<JSONOutputArchiveNoPrettyWriter>(this),
        itsWriteStream(stream),
        itsWriter(itsWriteStream),
        itsNextName(nullptr)
      {
        itsWriter.SetMaxDecimalPlaces( options.itsPrecision );
        itsNameCounter.push(0);
        itsNodeStack.push(NodeType::StartObject);
      }

      //! Destructor, flushes the JSON
      ~JSONOutputArchiveNoPrettyWriter() CEREAL_NOEXCEPT
      {
        if (itsNodeStack.top() == NodeType::InObject)
          itsWriter.EndObject();
        else if (itsNodeStack.top() == NodeType::InArray)
          itsWriter.EndArray();
      }

      //! Saves some binary data, encoded as a base64 string, with an optional name
      /*! This will create a new node, optionally named, and insert a value that consists of
          the data encoded as a base64 string */
      void saveBinaryValue( const void * data, size_t size, const char * name = nullptr )
      {
        setNextName( name );
        writeName();

        auto base64string = base64::encode( reinterpret_cast<const unsigned char *>( data ), size );
        saveValue( base64string );
      };

      //! @}
      /*! @name Internal Functionality
          Functionality designed for use by those requiring control over the inner mechanisms of
          the JSONOutputArchiveNoPrettyWriter */
      //! @{

      //! Starts a new node in the JSON output
      /*! The node can optionally be given a name by calling setNextName prior
          to creating the node

          Nodes only need to be started for types that are themselves objects or arrays */
      void startNode()
      {
        writeName();
        itsNodeStack.push(NodeType::StartObject);
        itsNameCounter.push(0);
      }

      //! Designates the most recently added node as finished
      void finishNode()
      {
        // if we ended up serializing an empty object or array, writeName
        // will never have been called - so start and then immediately end
        // the object/array.
        //
        // We'll also end any object/arrays we happen to be in
        switch(itsNodeStack.top())
        {
          case NodeType::StartArray:
            itsWriter.StartArray();
          case NodeType::InArray:
            itsWriter.EndArray();
            break;
          case NodeType::StartObject:
            itsWriter.StartObject();
          case NodeType::InObject:
            itsWriter.EndObject();
            break;
        }

        itsNodeStack.pop();
        itsNameCounter.pop();
      }

      //! Sets the name for the next node created with startNode
      void setNextName( const char * name )
      {
        itsNextName = name;
      }

      //! Saves a bool to the current node
      void saveValue(bool b)                { itsWriter.Bool(b);                                                         }
      //! Saves an int to the current node
      void saveValue(int i)                 { itsWriter.Int(i);                                                          }
      //! Saves a uint to the current node
      void saveValue(unsigned u)            { itsWriter.Uint(u);                                                         }
      //! Saves an int64 to the current node
      void saveValue(int64_t i64)           { itsWriter.Int64(i64);                                                      }
      //! Saves a uint64 to the current node
      void saveValue(uint64_t u64)          { itsWriter.Uint64(u64);                                                     }
      //! Saves a double to the current node
      void saveValue(double d)              { itsWriter.Double(d);                                                       }
      //! Saves a string to the current node
      void saveValue(std::string const & s) { itsWriter.String(s.c_str(), static_cast<CEREAL_RAPIDJSON_NAMESPACE::SizeType>( s.size() )); }
      //! Saves a const char * to the current node
      void saveValue(char const * s)        { itsWriter.String(s);                                                       }
      //! Saves a nullptr to the current node
      void saveValue(std::nullptr_t)        { itsWriter.Null();                                                          }

    private:
      // Some compilers/OS have difficulty disambiguating the above for various flavors of longs, so we provide
      // special overloads to handle these cases.

      //! 32 bit signed long saving to current node
      template <class T, traits::EnableIf<sizeof(T) == sizeof(std::int32_t),
                                          std::is_signed<T>::value> = traits::sfinae> inline
      void saveLong(T l){ saveValue( static_cast<std::int32_t>( l ) ); }

      //! non 32 bit signed long saving to current node
      template <class T, traits::EnableIf<sizeof(T) != sizeof(std::int32_t),
                                          std::is_signed<T>::value> = traits::sfinae> inline
      void saveLong(T l){ saveValue( static_cast<std::int64_t>( l ) ); }

      //! 32 bit unsigned long saving to current node
      template <class T, traits::EnableIf<sizeof(T) == sizeof(std::int32_t),
                                          std::is_unsigned<T>::value> = traits::sfinae> inline
      void saveLong(T lu){ saveValue( static_cast<std::uint32_t>( lu ) ); }

      //! non 32 bit unsigned long saving to current node
      template <class T, traits::EnableIf<sizeof(T) != sizeof(std::int32_t),
                                          std::is_unsigned<T>::value> = traits::sfinae> inline
      void saveLong(T lu){ saveValue( static_cast<std::uint64_t>( lu ) ); }

    public:
#ifdef _MSC_VER
      //! MSVC only long overload to current node
      void saveValue( unsigned long lu ){ saveLong( lu ); };
#else // _MSC_VER
      //! Serialize a long if it would not be caught otherwise
      template <class T, traits::EnableIf<std::is_same<T, long>::value,
                                          !std::is_same<T, std::int32_t>::value,
                                          !std::is_same<T, std::int64_t>::value> = traits::sfinae> inline
      void saveValue( T t ){ saveLong( t ); }

      //! Serialize an unsigned long if it would not be caught otherwise
      template <class T, traits::EnableIf<std::is_same<T, unsigned long>::value,
                                          !std::is_same<T, std::uint32_t>::value,
                                          !std::is_same<T, std::uint64_t>::value> = traits::sfinae> inline
      void saveValue( T t ){ saveLong( t ); }
#endif // _MSC_VER

      //! Save exotic arithmetic as strings to current node
      /*! Handles long long (if distinct from other types), unsigned long (if distinct), and long double */
      template <class T, traits::EnableIf<std::is_arithmetic<T>::value,
                                          !std::is_same<T, long>::value,
                                          !std::is_same<T, unsigned long>::value,
                                          !std::is_same<T, std::int64_t>::value,
                                          !std::is_same<T, std::uint64_t>::value,
                                          (sizeof(T) >= sizeof(long double) || sizeof(T) >= sizeof(long long))> = traits::sfinae> inline
      void saveValue(T const & t)
      {
        std::stringstream ss; ss.precision( std::numeric_limits<long double>::max_digits10 );
        ss << t;
        saveValue( ss.str() );
      }

      //! Write the name of the upcoming node and prepare object/array state
      /*! Since writeName is called for every value that is output, regardless of
          whether it has a name or not, it is the place where we will do a deferred
          check of our node state and decide whether we are in an array or an object.

          The general workflow of saving to the JSON archive is:

            1. (optional) Set the name for the next node to be created, usually done by an NVP
            2. Start the node
            3. (if there is data to save) Write the name of the node (this function)
            4. (if there is data to save) Save the data (with saveValue)
            5. Finish the node
          */
      void writeName()
      {
        NodeType const & nodeType = itsNodeStack.top();

        // Start up either an object or an array, depending on state
        if(nodeType == NodeType::StartArray)
        {
          itsWriter.StartArray();
          itsNodeStack.top() = NodeType::InArray;
        }
        else if(nodeType == NodeType::StartObject)
        {
          itsNodeStack.top() = NodeType::InObject;
          itsWriter.StartObject();
        }

        // Array types do not output names
        if(nodeType == NodeType::InArray) return;

        if(itsNextName == nullptr)
        {
          std::string name = "value" + std::to_string( itsNameCounter.top()++ ) + "\0";
          saveValue(name);
        }
        else
        {
          saveValue(itsNextName);
          itsNextName = nullptr;
        }
      }

      //! Designates that the current node should be output as an array, not an object
      void makeArray()
      {
        itsNodeStack.top() = NodeType::StartArray;
      }

      //! @}

    private:
      WriteStream itsWriteStream;          //!< Rapidjson write stream
      JSONWriter itsWriter;                //!< Rapidjson writer
      char const * itsNextName;            //!< The next name
      std::stack<uint32_t> itsNameCounter; //!< Counter for creating unique names for unnamed nodes
      std::stack<NodeType> itsNodeStack;
  }; // JSONOutputArchiveNoPrettyWriter

  // ######################################################################
  // JSONArchive prologue and epilogue functions
  // ######################################################################

  // ######################################################################
  //! Prologue for NVPs for JSON archives
  /*! NVPs do not start or finish nodes - they just set up the names */
  template <class T> inline
  void prologue( JSONOutputArchiveNoPrettyWriter &, NameValuePair<T> const & )
  { }

  // ######################################################################
  //! Epilogue for NVPs for JSON archives
  /*! NVPs do not start or finish nodes - they just set up the names */
  template <class T> inline
  void epilogue( JSONOutputArchiveNoPrettyWriter &, NameValuePair<T> const & )
  { }

  // ######################################################################
  //! Prologue for SizeTags for JSON archives
  /*! SizeTags are strictly ignored for JSON, they just indicate
      that the current node should be made into an array */
  template <class T> inline
  void prologue( JSONOutputArchiveNoPrettyWriter & ar, SizeTag<T> const & )
  {
    ar.makeArray();
  }

  // ######################################################################
  //! Epilogue for SizeTags for JSON archives
  /*! SizeTags are strictly ignored for JSON */
  template <class T> inline
  void epilogue( JSONOutputArchiveNoPrettyWriter &, SizeTag<T> const & )
  { }

  // ######################################################################
  //! Prologue for all other types for JSON archives (except minimal types)
  /*! Starts a new node, named either automatically or by some NVP,
      that may be given data by the type about to be archived

      Minimal types do not start or finish nodes */
  template <class T, traits::EnableIf<!std::is_arithmetic<T>::value,
                                      !traits::has_minimal_base_class_serialization<T, traits::has_minimal_output_serialization, JSONOutputArchiveNoPrettyWriter>::value,
                                      !traits::has_minimal_output_serialization<T, JSONOutputArchiveNoPrettyWriter>::value> = traits::sfinae>
  inline void prologue( JSONOutputArchiveNoPrettyWriter & ar, T const & )
  {
    ar.startNode();
  }

  // ######################################################################
  //! Epilogue for all other types other for JSON archives (except minimal types)
  /*! Finishes the node created in the prologue

      Minimal types do not start or finish nodes */
  template <class T, traits::EnableIf<!std::is_arithmetic<T>::value,
                                      !traits::has_minimal_base_class_serialization<T, traits::has_minimal_output_serialization, JSONOutputArchiveNoPrettyWriter>::value,
                                      !traits::has_minimal_output_serialization<T, JSONOutputArchiveNoPrettyWriter>::value> = traits::sfinae>
  inline void epilogue( JSONOutputArchiveNoPrettyWriter & ar, T const & )
  {
    ar.finishNode();
  }

  // ######################################################################
  //! Prologue for arithmetic types for JSON archives
  inline
  void prologue( JSONOutputArchiveNoPrettyWriter & ar, std::nullptr_t const & )
  {
    ar.writeName();
  }

  // ######################################################################
  //! Epilogue for arithmetic types for JSON archives
  inline
  void epilogue( JSONOutputArchiveNoPrettyWriter &, std::nullptr_t const & )
  { }

  // ######################################################################
  //! Prologue for arithmetic types for JSON archives
  template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae> inline
  void prologue( JSONOutputArchiveNoPrettyWriter & ar, T const & )
  {
    ar.writeName();
  }

  // ######################################################################
  //! Epilogue for arithmetic types for JSON archives
  template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae> inline
  void epilogue( JSONOutputArchiveNoPrettyWriter &, T const & )
  { }

  // ######################################################################
  //! Prologue for strings for JSON archives
  template<class CharT, class Traits, class Alloc> inline
  void prologue(JSONOutputArchiveNoPrettyWriter & ar, std::basic_string<CharT, Traits, Alloc> const &)
  {
    ar.writeName();
  }

  // ######################################################################
  //! Epilogue for strings for JSON archives
  template<class CharT, class Traits, class Alloc> inline
  void epilogue(JSONOutputArchiveNoPrettyWriter &, std::basic_string<CharT, Traits, Alloc> const &)
  { }

  // ######################################################################
  // Common JSONArchive serialization functions
  // ######################################################################
  //! Serializing NVP types to JSON
  template <class T> inline
  void CEREAL_SAVE_FUNCTION_NAME( JSONOutputArchiveNoPrettyWriter & ar, NameValuePair<T> const & t )
  {
    ar.setNextName( t.name );
    ar( t.value );
  }

  //! Saving for nullptr to JSON
  inline
  void CEREAL_SAVE_FUNCTION_NAME(JSONOutputArchiveNoPrettyWriter & ar, std::nullptr_t const & t)
  {
    ar.saveValue( t );
  }

  //! Saving for arithmetic to JSON
  template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae> inline
  void CEREAL_SAVE_FUNCTION_NAME(JSONOutputArchiveNoPrettyWriter & ar, T const & t)
  {
    ar.saveValue( t );
  }

  //! saving string to JSON
  template<class CharT, class Traits, class Alloc> inline
  void CEREAL_SAVE_FUNCTION_NAME(JSONOutputArchiveNoPrettyWriter & ar, std::basic_string<CharT, Traits, Alloc> const & str)
  {
    ar.saveValue( str );
  }

  // ######################################################################
  //! Saving SizeTags to JSON
  template <class T> inline
  void CEREAL_SAVE_FUNCTION_NAME( JSONOutputArchiveNoPrettyWriter &, SizeTag<T> const & )
  {
    // nothing to do here, we don't explicitly save the size
  }

} // namespace cereal

// register archives for polymorphic support
CEREAL_REGISTER_ARCHIVE(cereal::JSONOutputArchiveNoPrettyWriter)

#endif // CEREAL_ARCHIVES_JSON_NO_PRETTY_WRITER_HPP_
