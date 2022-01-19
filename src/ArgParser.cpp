#include <map>
#include <string>
#include <iostream>

typedef std::map<std::string, std::string>::const_iterator it_type;

class ArgParser
{

	public:

	ArgParser(const int argc, char* const argv[])
  {
		for (int i = 0; i < argc; i++)
		{
			std::string str = argv[i]; 
			if (str[0] != '-') 
				continue;
			if (str[1] == '-')
				index[str.substr(2, str.length() - 2)] = " ";
			else if (i < argc - 1)	
				index[str.substr(1, str.length() - 1)] = argv[++i];  		
		}
  }

  std::string operator()(const std::string key)
  {
    if (exists(key)) 
    	return index[key];
    return "";	
  }
  
  bool exists(const std::string key)
  {
  	return index.count(key);
  }

	private:

  std::map<std::string, std::string> index;
	
	friend std::ostream& operator<<(std::ostream& os, const ArgParser& ap)
  {
		for(auto& it : ap.index) 
		{
			os << it.first << " " << it.second << std::endl;
		}
    return os;
  }

  ArgParser() { }
  
};
