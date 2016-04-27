// Header Include
#include "CommandLineParser.h"
#include "Logger.h"

// System Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

namespace specnets
{

  //----------------------------------------------------------------------------
  CommandLineParser::CommandLineParser(int argc,
                                       char ** argv,
                                       int numNonOptions,
                                       const std::vector<Option> & listOptions)
  {
    for (size_t i = numNonOptions + 1; i < argc; i++)
    {
      //DEBUG_VAR(i);
      //DEBUG_VAR(argv[i]);

      bool unknownOption = true;
      //DEBUG_VAR(listOptions.size());
      for (size_t j = 0; j < listOptions.size(); j++)
      {
        //DEBUG_VAR(j);

        string option("-");
        option += listOptions[j].m_flag;
        //DEBUG_VAR(option);

        if (::strcmp(argv[i], option.c_str()) == 0)
        {
          //DEBUG_TRACE;
          unknownOption = false;
          //DEBUG_VAR(listOptions[j].m_numValues);
          if (listOptions[j].m_numValues != 0)
          {
            //DEBUG_TRACE;
            for (unsigned int v = 0; v < listOptions[j].m_numValues; v++)
            {
              i++;
              if (i == argc)
              {
                //DEBUG_TRACE;
                char buf[256];
                sprintf(buf, "%u", listOptions[j].m_numValues);
                m_error = "Option [" + option + "] requires [" + buf + "] argument(s)";
                return;
              }
              else
              {
                //DEBUG_VAR(listOptions[j].m_name);
                //DEBUG_VAR(argv[i]);
                m_listOptions.push_back(make_pair<string, string> (listOptions[j].m_name,
                                                                   argv[i]));
              }
            }
          }
          else
          {
            //DEBUG_VAR(listOptions[j].m_name);
            m_listOptions.push_back(make_pair<string, string> (listOptions[j].m_name,
                                                               ""));
          }
        }

      } // for (size_t j = 0; j < listOptions.size(); j++)

      if (unknownOption)
      {
        m_error = "Unknown option [";
        m_error += argv[i];
        m_error += "]";
        return;
      }

    } // for (size_t i = numNonOptions + 1; i < argc; i++)

    return;
  }

  CommandLineParser::~CommandLineParser()
  {
    // EMPTY
  }

  //----------------------------------------------------------------------------
  void CommandLineParser::getOptionsAsParameterList(ParameterList & params) const
  {
    for (size_t i = 0; i < m_listOptions.size(); i++)
    {
      //DEBUG_VAR(m_listOptions[i].first);
      //DEBUG_VAR(m_listOptions[i].second);
      params.setValue(m_listOptions[i].first, m_listOptions[i].second);
    }
    return;
  }

  //----------------------------------------------------------------------------
  bool CommandLineParser::validate(string & error) const
  {
    if (m_error.empty())
    {
      return true;
    }
    error = m_error;
    return false;
  }

} //namespace specnets
