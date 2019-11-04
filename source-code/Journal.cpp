#ifndef C_JOURNAL
#define C_JOURNAL

#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <iomanip>

struct item {

  std::string item_s;
  struct item *next_item;
  struct item *prev_item;

};

class Journal {

  private:
    std::ofstream f_journal;
    int offset;
    char *ReadTime();
    struct item *first_item;
    struct item *last_item;
    int count_items;

  protected:


  public:

    Journal();
    Journal(std::string in_file);
    ~Journal();

    int JournalIn(std::string in_func);
    int JournalW(std::string message);
    int JournalOut();

};

Journal :: Journal() {

  this->offset = 0;
  this->first_item = this->last_item = NULL;
  count_items = 0;
  std::cout << "Journal init" << std::endl;

}

Journal :: Journal(std::string in_file) {

  this->f_journal.open(in_file.c_str());
  this->offset = 0;
  this->first_item = this->last_item = NULL;
  count_items = 0;
  this->f_journal << "Journal opened in file " << in_file << " at " << this->ReadTime() << std::endl;

}

Journal :: ~Journal() {

  this->f_journal << "Journal closed at " << this->ReadTime() << std::endl;
  this->f_journal.close();

}

int Journal :: JournalW(std::string message) {

  int i;
  for(i = 0; i < this->offset / 2; i++) this->f_journal << " |";

  this->f_journal << " |" << message << std::setw(80 - message.length() - offset - 2) << this->ReadTime() << std::endl;

  return(0);

}

int Journal :: JournalIn(std::string in_func) {

  int i;
  this->offset += 2;
  for(i = 0; i < this->offset / 2; i++) this->f_journal << " |";
  this->f_journal << in_func << " start" << std::setw(80 - in_func.length() - offset - 6) << this->ReadTime() << std::endl;
  if(count_items == 0) {
    this->first_item = this->last_item = new item;
    first_item->next_item = first_item->prev_item = NULL;
    first_item->item_s = in_func;
  }
  else {
    last_item->next_item = new item;
    last_item->next_item->prev_item = last_item;
    last_item->next_item->next_item = NULL;
    last_item = last_item->next_item;
    last_item->item_s = in_func;
  }
  count_items++;
  return(0);

}

int Journal :: JournalOut() {

  int i;
  for(i = 0; i < this->offset / 2; i++) this->f_journal << " |";
  this->f_journal << last_item->item_s << " end" << std::setw(80 - last_item->item_s.length() - offset - 4) << this->ReadTime() << std::endl;
  this->offset -= 2;
  if(count_items > 0) {
    if(count_items > 1) {
      last_item = last_item->prev_item;
      delete last_item->next_item;
      last_item->next_item = NULL;
    } else {
      delete last_item;
      first_item = last_item = NULL;
    }
    count_items--;
  }
  return(0);

}

char *Journal :: ReadTime() {

  time_t now = time(0);
  char *datetime = ctime(&now);
  datetime[24] = '\0';
  return(datetime);

}

#endif // C_JOURNAL
