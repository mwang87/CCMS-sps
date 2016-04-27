/*
 * DelimitedTextReader.cpp
 *
 *  Created on: Jan 2, 2011
 *      Author: jsnedecor
 */

#include "DelimitedTextReader.h"

namespace specnets
{
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::checkHeader(map<string, unsigned int> &header,
                                        vector<string> &requiredHeader,
                                        vector<int> &requiredHeaderIndex,
                                        int &missingIndex)
  {
    bool allMatched = true; //check to see whether we've matched required headers

    map<string, unsigned int>::const_iterator headerIter;

    requiredHeaderIndex.resize(requiredHeader.size());

    for (unsigned int i = 0; i < requiredHeader.size(); i++)
    {
      headerIter = header.find(requiredHeader[i]);
      if (headerIter == header.end())
      {
        missingIndex = i;
        allMatched = false;
      }
      else
      {
        requiredHeaderIndex[i] = headerIter->second;
      }
    }
    return allMatched;
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::checkHeader(vector<string> &header,
                                        vector<string> &requiredHeader,
                                        vector<int> &requiredHeaderIndex,
                                        int &missingIndex)
  {
    bool allMatched = true; //check to see whether we've matched required headers

    requiredHeaderIndex.resize(requiredHeader.size());

    for (unsigned int i = 0; i < requiredHeader.size(); i++)
    {
      bool currMatched = false;
      for (unsigned int j = 0; j < header.size(); j++)
      {
        if (header[j].compare(requiredHeader[i]) == 0)
        {
          currMatched = true;
          requiredHeaderIndex[i] = j;
        }
      }
      if (!currMatched)
      {
        missingIndex = i;
        allMatched = false;
      }
    }
    return allMatched;
  }
  //---------------------------------------------------------------------------
  template<typename T> bool loadDelimitedFileNoHeaderT(const char * filename,
                                                       const string &delim,
                                                       const string &commentDelim,
                                                       vector<vector<T> > &lines)
  {
    //open delimited file
    ifstream delimitedFileHandle(filename);

    if (!delimitedFileHandle.is_open() || !delimitedFileHandle.good())
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }

    //read fields for file
    string lineBuff;

    while (!delimitedFileHandle.eof())
    {
      DelimitedTextReader::getlineChomp(delimitedFileHandle, lineBuff);

      if (lineBuff[0] == *commentDelim.c_str())
      { //skip comments
        //do nothing
        DEBUG_TRACE;
      }
      else if (lineBuff.compare("") == 0) //skip blank lines
      {
        DEBUG_TRACE;
      }
      else
      {
        //load in  file
        vector<string> fields;
        stringSplit2(lineBuff, fields, delim);
        insertValue(fields, lines);
      }
    }
    return true;
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(const char * filename,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<float> > &lines)
  {
    return loadDelimitedFileNoHeaderT(filename, delim, commentDelim, lines);
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(const char * filename,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<int> > &lines)
  {
    return loadDelimitedFileNoHeaderT(filename, delim, commentDelim, lines);
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(const char * filename,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<string> > &lines)
  {
    return loadDelimitedFileNoHeaderT(filename, delim, commentDelim, lines);
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(const char * filename,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<double> > &lines)
  {
    return loadDelimitedFileNoHeaderT(filename, delim, commentDelim, lines);
  }

  //---------------------------------------------------------------------------
  template<typename T> bool loadDelimitedFileNoHeaderT(ifstream &delimitedFileHandle,
                                                       const string &delim,
                                                       const string &commentDelim,
                                                       vector<vector<T> > &lines)
  {
    //read fields for file
    string lineBuff;

    while (!delimitedFileHandle.eof())
    {
      DelimitedTextReader::getlineChomp(delimitedFileHandle, lineBuff);

      if (lineBuff[0] == *commentDelim.c_str())
      { //skip comments
        //do nothing
        DEBUG_TRACE;
      }
      else
      {
        //load in  file
        vector<string> fields;
        stringSplit2(lineBuff, fields, delim);
        insertValue(fields, lines);
      }
    }
    return true;
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(ifstream &delimitedFileHandle,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<double> > &lines)
  {
    return loadDelimitedFileNoHeaderT(delimitedFileHandle,
                                      delim,
                                      commentDelim,
                                      lines);
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(ifstream &delimitedFileHandle,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<float> > &lines)
  {
    return loadDelimitedFileNoHeaderT(delimitedFileHandle,
                                      delim,
                                      commentDelim,
                                      lines);
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(ifstream &delimitedFileHandle,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<string> > &lines)
  {
    return loadDelimitedFileNoHeaderT(delimitedFileHandle,
                                      delim,
                                      commentDelim,
                                      lines);
  }

  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFileNoHeader(ifstream &delimitedFileHandle,
                                                      const string &delim,
                                                      const string &commentDelim,
                                                      vector<vector<int> > &lines)
  {
    return loadDelimitedFileNoHeaderT(delimitedFileHandle,
                                      delim,
                                      commentDelim,
                                      lines);
  }
  //---------------------------------------------------------------------------

  template<typename T> bool loadDelimitedFileT(const char * filename,
                                               const string &delim,
                                               const string &commentDelim,
                                               map<string, unsigned int> &header,
                                               vector<vector<T> > &lines,
                                               vector<string> &requiredHeader,
                                               vector<int> &requiredHeaderIndex)
  {
    //open delimited file
    ifstream delimitedFileHandle(filename, ios::binary);

    if (!delimitedFileHandle.is_open() || !delimitedFileHandle.good())
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }

    //read fields for file
    string lineBuff;

    bool readHeader = false;

    while (!delimitedFileHandle.eof())
    {
      DelimitedTextReader::getlineChomp(delimitedFileHandle, lineBuff);

      //skip blank lines
      if (lineBuff.compare("") == 0)
      {
        continue;
      }

      if (lineBuff[0] == *commentDelim.c_str())
      { //skip comments
        //do nothing
        DEBUG_TRACE;
      }
      else
      {
        if (!readHeader)
        { //haven't yet gotten header
          vector<string> headerTemp;
          stringSplit2(lineBuff, headerTemp, delim);

          for (unsigned int i = 0; i < headerTemp.size(); i++)
          {
            header[headerTemp[i]] = i;
            // DEBUG_MSG("header " << headerTemp[i] << delim << i);
          }

          readHeader = true;

          int missingIndex;

          if (!DelimitedTextReader::checkHeader(header,
                                                requiredHeader,
                                                requiredHeaderIndex,
                                                missingIndex))
          {
            ERROR_MSG("Missing required header! " << requiredHeader[missingIndex]);
            return false;
          }
          // DEBUG_MSG("Checked header");
        }
        else
        {
          //load in  file, header already defined
          vector<string> fields;
          stringSplit2(lineBuff, fields, delim);
          insertValue(fields, lines);
        }
      }
    }
    return true;
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              map<string, unsigned int> &header,
                                              vector<vector<int> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              map<string, unsigned int> &header,
                                              vector<vector<float> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              map<string, unsigned int> &header,
                                              vector<vector<double> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              map<string, unsigned int> &header,
                                              vector<vector<string> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }

  //---------------------------------------------------------------------------

  template<typename T> bool loadDelimitedFileT(const char * filename,
                                               const string &delim,
                                               const string &commentDelim,
                                               vector<string> &header,
                                               vector<vector<T> > &lines,
                                               vector<string> &requiredHeader,
                                               vector<int> &requiredHeaderIndex)
  {
    //open delimited file
    ifstream delimitedFileHandle(filename);

    if (!delimitedFileHandle.is_open() || !delimitedFileHandle.good())
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }

    //read fields for file
    string lineBuff;

    bool readHeader = false;

    while (!delimitedFileHandle.eof() && delimitedFileHandle.good())
    {
      DelimitedTextReader::getlineChomp(delimitedFileHandle, lineBuff);

      //skip blank lines
      if (lineBuff.compare("") == 0)
      {
        continue;
      }

      if (lineBuff[0] == *commentDelim.c_str())
      { //skip comments
        //do nothing
        DEBUG_TRACE;
      }
      else
      {
        if (!readHeader)
        { //haven't yet gotten header
          stringSplit2(lineBuff, header, delim);

          readHeader = true;

          int missingIndex;

          if (!DelimitedTextReader::checkHeader(header,
                                                requiredHeader,
                                                requiredHeaderIndex,
                                                missingIndex))
          {
            ERROR_MSG("Missing required header! " << requiredHeader[missingIndex]);
            return false;
          }
          // DEBUG_MSG("Checked header");
        }
        else
        {
          //load in  file, header already defined
          vector<string> fields;
          stringSplit2(lineBuff, fields, delim);
          insertValue(fields, lines);
        }
      }
    }
    return true;
  }
  //---------------------------------------------------------------------------

  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              vector<string> &header,
                                              vector<vector<string> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              vector<string> &header,
                                              vector<vector<float> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              vector<string> &header,
                                              vector<vector<double> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  bool DelimitedTextReader::loadDelimitedFile(const char * filename,
                                              const string &delim,
                                              const string &commentDelim,
                                              vector<string> &header,
                                              vector<vector<int> > &lines,
                                              vector<string> &requiredHeader,
                                              vector<int> &requiredHeaderIndex)
  {
    return loadDelimitedFileT(filename,
                              delim,
                              commentDelim,
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
  }
  //---------------------------------------------------------------------------
  void DelimitedTextReader::getlineChomp(ifstream &fileHandle,
                                         string &lineBuffer)
  {
    if (getline(fileHandle, lineBuffer))
    {
      if ('\r' == lineBuffer[lineBuffer.size() - 1])
      { //strip off carriage return for dos
        lineBuffer.resize(lineBuffer.size() - 1);
      }
    }
  }
  //---------------------------------------------------------------------------

  void insertValue(vector<string> & fields, vector<vector<string> > &lines)
  {
    lines.push_back(fields);
  }

  //---------------------------------------------------------------------------
  void insertValue(vector<string> & fields, vector<vector<float> > &lines)
  {
    vector<float> temp;

    for (int i = 0; i < fields.size(); i++)
    {
      float curr_field;
      sscanf(fields[i].c_str(), "%f", &curr_field);
      temp.push_back(curr_field);
    }
    lines.push_back(temp);
  }

  //---------------------------------------------------------------------------
  void insertValue(vector<string> & fields, vector<vector<int> > &lines)
  {
    vector<int> temp;

    for (int i = 0; i < fields.size(); i++)
    {
      int curr_field;
      sscanf(fields[i].c_str(), "%d", &curr_field);
      temp.push_back(curr_field);
    }
    lines.push_back(temp);
  }

  //---------------------------------------------------------------------------
  void insertValue(vector<string> & fields, vector<vector<double> > &lines)
  {
    vector<double> temp;

    for (int i = 0; i < fields.size(); i++)
    {
      double curr_field;
      sscanf(fields[i].c_str(), "%lf", &curr_field);
      temp.push_back(curr_field);
    }

    lines.push_back(temp);
  }
}
