import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import check_phase_a_reuse as reuse


class PhaseAReuseTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.jobs=[]
        for job_id in reuse.EXPECTED_JOB_ORDER[:-1]:
            self.jobs.append({"id":job_id,"status":"completed","returncode":0,
                "inputs_unchanged":True,"cleanup":{"verified_gone":True},"failure_reasons":[]})
        log=self.root/reuse.A4_GATE_ID/"job.log"
        log.parent.mkdir();log.write_text("ValueError: A4 profiler call/RHS counters failed\n")
        self.jobs.append({"id":reuse.A4_GATE_ID,"status":"failed","returncode":1,
            "inputs_unchanged":True,"cleanup":{"verified_gone":True},
            "failure_reasons":["missing_outputs:checks.json"],"log":"job.log",
            "log_sha256":hashlib.sha256(log.read_bytes()).hexdigest()})
        self.state_path=self.root/"state.json"
        self.write_state()

    def tearDown(self):
        self.temp.cleanup()

    def write_state(self):
        self.state_path.write_text(json.dumps({"status":"failed","jobs":self.jobs},separators=(",",":")))
        self.state_hash=hashlib.sha256(self.state_path.read_bytes()).hexdigest()

    def test_requalifies_only_the_known_a4_checker_failure(self):
        with patch.object(reuse,"TASK_STATE_SHA256",self.state_hash), \
             patch.object(reuse.check_pair,"verified_job",return_value=({"returncode":0},{"data/a4-summary.json":"a"},self.state_hash)), \
             patch.object(reuse.check_phase_a,"data_gate",return_value={"status":"pass","relative_l2":0.0}):
            result=reuse.qualify(self.root,self.state_hash)
        self.assertEqual(result["status"],"pass")
        self.assertEqual(result["rechecked_a4_gate"]["status"],"pass")

    def test_rejects_any_a5_job_in_the_failed_source_task(self):
        self.jobs.append({"id":"A5-profiled-standard-forward-two-submodes"})
        self.write_state()
        with patch.object(reuse,"TASK_STATE_SHA256",self.state_hash), \
             self.assertRaisesRegex(ValueError,"inventory/order changed or A5 already started"):
            reuse.qualify(self.root,self.state_hash)


if __name__=="__main__":
    unittest.main()
