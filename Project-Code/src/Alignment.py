
class Alignment:
  
  def __init__(self,hit_id, query_range, hit_range, query_seq, hit_seq, midline):
    self.hit_id = hit_id
    self.query_range = query_range 
    self.hit_range = hit_range

    self.query_seq = query_seq
    self.hit_seq = hit_seq
    self.midline = midline

    self.query_hit_dict = self.generate_query_hit_dict()


  def generate_query_hit_dict(self):
    query_hit_dict = {}
    q_index, h_index = self.query_range[0], self.hit_range[0]
    for q, h, m in zip(self.query_seq, self.hit_seq, self.midline):
      # if gap in query sequence
      if q == '-':
        h_index += 1
      # if gap in hit sequence
      elif h == '-':
        q_index += 1
      # if residues align
      elif m != ' ':
        query_hit_dict[q_index] = h_index
        q_index += 1
        h_index += 1
      # if residues do not align
      else:
        q_index += 1
        h_index += 1

    return query_hit_dict